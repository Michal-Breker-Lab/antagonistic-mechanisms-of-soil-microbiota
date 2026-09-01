#!/usr/bin/env bash
# Add a SECOND annotation worker on the far end of the job list.
#
# Why not just restart stage 2 with a higher -P: that means killing 14 in-flight
# genomes, and pkill -f over ssh matches the ssh command line and kills the
# calling session. Nothing is killed here.
#
# Race avoidance: worker 1 walks the list top-down at ~24 genomes/h. This worker
# takes only the BOTTOM half, so worker 1 cannot reach this worker's territory
# for ~9 h, by which time these genomes carry a .gff3 and worker 1's run_one
# skips them. Convergence is monitored; this worker is stopped by PIDFILE well
# before the two meet.
set -uo pipefail
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
export TMPDIR=$W/tmp
DB=/mnt/LargeStorageNoBackup/Datasets/Moshea/Databases/bakta_db/db
BAKTA=$W/envs/bakta/bin/bakta
export PATH=$W/envs/bakta/bin:$PATH
JOBS=${JOBS:-14}
THREADS=${THREADS:-8}
SRC=$W/annot_jobs_full.tsv
MINE=$W/annot_jobs_worker2.tsv

[ -s "$SRC" ] || { echo "FATAL: $SRC missing" >&2; exit 1; }
TOTAL=$(wc -l < "$SRC")
HALF=$(( TOTAL / 2 ))
# bottom half, minus anything already finished
tail -n "$HALF" "$SRC" | while IFS=$'\t' read -r acc fna genus sp comp; do
  [ -s "$W/annot/$acc/$acc.gff3" ] || printf '%s\t%s\t%s\t%s\t%s\n' "$acc" "$fna" "$genus" "$sp" "$comp"
done > "$MINE"

echo "job list total     : $TOTAL"
echo "worker 2 territory : bottom $HALF"
echo "still to do here   : $(wc -l < "$MINE")"
echo "first in territory : $(head -1 "$MINE" | cut -f1)"
echo "concurrency        : $JOBS x $THREADS threads"
[ "${1:-}" = "--confirm" ] || { echo; echo "Dry run. Re-run with --confirm."; exit 0; }

echo $$ > "$W/logs/worker2.pid"
run_one() {
  IFS=$'\t' read -r acc fna genus sp comp <<<"$1"
  out=$W/annot/$acc
  # re-check under the wire: worker 1 may have taken it since the list was built
  if [ -s "$out/$acc.gff3" ]; then echo "skip $acc"; return 0; fi
  mkdir -p "$out"
  cflag=""; [ "$comp" = "1" ] && cflag="--complete"
  $BAKTA --db $DB --threads $THREADS --genus "$genus" --species "$sp" \
         --prefix "$acc" --output "$out" $cflag --keep-contig-headers \
         --skip-plot --force --tmp-dir $TMPDIR "$fna" \
         > $W/logs/bakta/$acc.log 2>&1
  if [ -s "$out/$acc.gff3" ]; then echo "ok $acc"; else echo "FAIL $acc"; fi
}
export -f run_one; export W BAKTA DB THREADS TMPDIR PATH
date
cat "$MINE" | xargs -d '\n' -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}
date
echo "WORKER2_DONE"
