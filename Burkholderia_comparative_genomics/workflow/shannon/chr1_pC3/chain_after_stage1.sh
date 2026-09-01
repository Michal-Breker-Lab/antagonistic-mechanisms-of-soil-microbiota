#!/usr/bin/env bash
# Wait for stage 1, then fetch any genome whose accession identity changed, then
# stop. The heavy stages (2 annotate, 5 pangenome, 6 trees) are NOT chained --
# stage 2 needs its Bakta invocation reviewed and the machine needs free cores.
set -uo pipefail
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
LOG="$W/logs/chain1.log"
echo $$ > "$W/logs/chain1.pid"
exec >> "$LOG" 2>&1
echo "[$(date +%F' '%T)] chain waiting for stage 1"

for i in $(seq 1 720); do            # up to 2 h
  [ -f "$W/state/stage1.done" ] && break
  if grep -q "FAIL" "$W/logs/stage1.out" 2>/dev/null; then
    echo "[$(date +%F' '%T)] stage 1 FAILED - chain aborting"; exit 1
  fi
  sleep 10
done
if [ ! -f "$W/state/stage1.done" ]; then
  echo "[$(date +%F' '%T)] timed out waiting for stage 1"; exit 1
fi

echo "[$(date +%F' '%T)] stage 1 done; downloading any missing FASTA"
bash "$W/s1c_download_missing.sh"
rc=$?
echo "[$(date +%F' '%T)] s1c rc=$rc"
[ $rc -eq 0 ] || { echo "CHAIN_FAILED"; exit 1; }

echo "[$(date +%F' '%T)] genome set complete:"
echo "  list    : $(wc -l < "$W/genome_list_full.txt")"
echo "  fasta   : $(ls "$W"/genomes/*.fna 2>/dev/null | wc -l) + MF6 + MF7"
echo "  annot   : $(ls "$W"/annot 2>/dev/null | wc -l)"
touch "$W/state/chain1.done"
echo "CHAIN_OK"
