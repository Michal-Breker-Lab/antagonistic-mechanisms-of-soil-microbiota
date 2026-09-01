#!/usr/bin/env bash
# Read-only Shannon status probe. Safe to run any time; changes nothing.
# Reports whether the machine is free enough to start the rebuild.
set -uo pipefail
H=${H:-shannon}
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}

# NOTE: never use pgrep/pkill -f over ssh here - the pattern rides in the ssh
# command line, so pkill kills this very session (exit 255) and pgrep reports
# dead processes alive. Match with `ps -eo cmd | grep -E "[b]racket"` instead.
timeout 90 ssh -o BatchMode=yes -o ConnectTimeout=20 "$H" bash -s <<'REMOTE'
set -uo pipefail
cores=$(nproc)
read -r l1 l5 l15 rest < /proc/loadavg
printf "host        : %s\n" "$(hostname)"
printf "cores       : %s\n" "$cores"
printf "loadavg     : %s %s %s\n" "$l1" "$l5" "$l15"

# current idle from /proc/stat deltas -- NOT `ps %cpu`, which is a lifetime
# average and reports long-lived idle processes as busy.
read -r _ u1 n1 s1 i1 w1 rest1 < /proc/stat
sleep 2
read -r _ u2 n2 s2 i2 w2 rest2 < /proc/stat
dt=$(( (u2-u1)+(n2-n1)+(s2-s1)+(i2-i1)+(w2-w1) ))
di=$(( (i2-i1)+(w2-w1) ))
free_cores=$(awk -v di="$di" -v dt="$dt" -v c="$cores" 'BEGIN{printf "%.0f", (dt>0? di/dt*c : 0)}')
printf "free cores  : %s / %s\n" "$free_cores" "$cores"

printf "disk        : %s\n" "$(df -h /mnt/LargeStorageNoBackup | tail -1 | awk '{print $4" free ("$5" used)"}')"
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
printf "workdir     : %s\n" "$( [ -d "$W" ] && echo present || echo MISSING )"
printf "genomes     : %s\n" "$(ls "$W"/genomes 2>/dev/null | wc -l)"
printf "annotated   : %s\n" "$(ls "$W"/annot 2>/dev/null | wc -l)"
printf "MF7 staged  : %s\n" "$( [ -s "$W/MF7.fna" ] && echo yes || echo NO )"
printf "rebuild dir : %s\n" "$( [ -f "$W/run_full_rebuild.sh" ] && echo staged || echo NOT-STAGED )"

if [ "$free_cores" -ge 60 ]; then
  echo "VERDICT     : FREE - ok to launch"
else
  echo "VERDICT     : BUSY - hold"
  echo "top consumers:"
  ps -eo pcpu,rss,user,comm --sort=-pcpu | head -6 | sed 's/^/  /'
fi
REMOTE
rc=$?
[ $rc -ne 0 ] && echo "ssh failed (rc=$rc) - Shannon may be too loaded to accept commands"
exit $rc
