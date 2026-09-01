#!/usr/bin/env bash
# Download only the genomes that the full set names but genomes/ does not hold.
#
# Why not just re-run s1b_download.sh: it re-splits accessions.txt into
# batch_000, batch_001, ... and skips any batch whose .done marker exists. If the
# accession list has shifted by even one entry the batch boundaries move, the old
# markers no longer describe the new batches, and it silently skips the wrong
# accessions. This script never touches those markers -- it works from an
# explicit missing-list and downloads into its own directory.
set -uo pipefail
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
export TMPDIR=$W/tmp
DATASETS=$W/envs/burk/bin/datasets
LIST=${LIST:-$W/genome_list_full.txt}
mkdir -p "$W/raw/missing" "$W/genomes" "$W/logs" "$TMPDIR"

MISS="$W/logs/missing_fasta.txt"
: > "$MISS"
while read -r acc; do
  [ -z "$acc" ] && continue
  case "$acc" in MF6|MF7) continue;; esac        # lab isolates, not from NCBI
  [ -s "$W/genomes/$acc.fna" ] || echo "$acc" >> "$MISS"
done < "$LIST"

n=$(wc -l < "$MISS")
echo "full set          : $(wc -l < "$LIST")"
echo "already on disk   : $(ls "$W"/genomes/*.fna 2>/dev/null | wc -l)"
echo "missing FASTA     : $n"
[ "$n" -eq 0 ] && { echo "nothing to download"; exit 0; }
sed 's/^/  /' "$MISS"

$DATASETS download genome accession --inputfile "$MISS" \
    --include genome --no-progressbar \
    --filename "$W/raw/missing/missing.zip" > "$W/logs/dl_missing.log" 2>&1
rc=$?
# A truncated stream exits 0 and leaves a short zip, so check the log too.
if grep -qi "INTERNAL_ERROR\|stream error" "$W/logs/dl_missing.log"; then
  echo "FAIL: NCBI stream error - the zip is likely truncated, re-run"; exit 1
fi
[ $rc -eq 0 ] && [ -s "$W/raw/missing/missing.zip" ] || { echo "FAIL: download rc=$rc"; exit 1; }

unzip -o -q "$W/raw/missing/missing.zip" -d "$W/raw/missing/x" || { echo "FAIL: unzip"; exit 1; }
got=0
for gd in "$W"/raw/missing/x/ncbi_dataset/data/GC*/; do
  [ -d "$gd" ] || continue
  acc=$(basename "$gd")
  f=$(ls "$gd"*.fna 2>/dev/null | head -1)
  [ -n "$f" ] || continue
  cp -n "$f" "$W/genomes/$acc.fna"
  got=$((got+1))
done
echo "flattened         : $got"

still=0
while read -r acc; do
  [ -s "$W/genomes/$acc.fna" ] || { echo "  STILL MISSING $acc"; still=$((still+1)); }
done < "$MISS"
echo "still missing     : $still"
[ "$still" -eq 0 ] || exit 1
echo "OK"
