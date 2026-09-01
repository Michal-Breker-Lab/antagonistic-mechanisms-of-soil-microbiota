#!/usr/bin/env bash
# Smash-and-grab: the moment Shannon answers, pull the irreplaceable things FIRST.
#
# The s1-s27 pipeline scripts exist ONLY on Shannon. They are a few hundred KB
# and cannot be reconstructed without rewriting them, so they outrank everything
# else -- including the ~561 Bakta annotations, which are merely ~8 h of compute
# and can be regenerated on Moriah for about $12.
#
# Priority order matters because the window may be brief: Shannon has been
# oscillating between wedged and down for two days.
#
# tar over ssh, not rsync: rsync dies on every file over the DrvFs Drive mount.
set -uo pipefail
H=${H:-shannon}
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
D="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$D/from_shannon"
mkdir -p "$OUT"/{scripts,metadata,results,logs}

try() { timeout "${2:-120}" ssh -o BatchMode=yes -o ConnectTimeout=30 "$H" "$1" 2>/dev/null; }

for i in $(seq 1 500); do
  if try "true" 20; then
    echo "=== REACHED $H at $(date +%F' '%T) — grabbing ==="

    # 1. SCRIPTS + QUERY SEQUENCES (tiny, irreplaceable)
    echo "[1/4] pipeline scripts + queries"
    try "cd $W && tar cf - --ignore-failed-read \
         *.py *.sh *.R *.faa *.tsv *.txt 2>/dev/null" 300 \
      | tar xf - -C "$OUT/scripts" 2>/dev/null
    echo "      got $(ls -1 "$OUT/scripts" 2>/dev/null | wc -l) files"

    # 2. METADATA (few MB) — assemblies.tsv drives genus/species for annotation
    echo "[2/4] metadata"
    try "cd $W && tar cf - metadata metadata_20260803 2>/dev/null" 300 \
      | tar xf - -C "$OUT/metadata" 2>/dev/null

    # 3. RESULTS TABLES (tens of MB)
    echo "[3/4] results tables"
    try "cd $W/results && tar cf - --ignore-failed-read *.tsv *.txt 2>/dev/null" 600 \
      | tar xf - -C "$OUT/results" 2>/dev/null

    # 4. STATE: how far did annotation actually get?
    echo "[4/4] annotation state"
    try "ls -1 $W/annot/*/*.gff3 2>/dev/null | wc -l; uptime -p; \
         grep -c '^ok ' $W/logs/stage2.out 2>/dev/null; \
         grep -c '^FAIL' $W/logs/stage2.out 2>/dev/null" 120 \
      > "$OUT/logs/annotation_state.txt"
    echo "      $(cat "$OUT/logs/annotation_state.txt" | tr '\n' ' ')"

    echo "=== GRAB COMPLETE ==="
    du -sh "$OUT"
    exit 0
  fi
  sleep 60
done
echo "never reached $H in 500 attempts"
exit 1
