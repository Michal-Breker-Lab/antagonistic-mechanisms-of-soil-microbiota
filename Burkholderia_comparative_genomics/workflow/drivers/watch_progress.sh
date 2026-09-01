#!/usr/bin/env bash
# Poll annotation progress until stage 2 finishes (or the horizon expires).
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
for i in $(seq 1 200); do
  out=$(timeout 45 ssh -o BatchMode=yes -o ConnectTimeout=20 shannon "
      done_n=\$(ls $W/annot/*/*.gff3 2>/dev/null | wc -l)
      procs=\$(ps -eo cmd | grep -cE '[b]akta')
      fails=\$(grep -c '^FAIL' $W/logs/stage2.out 2>/dev/null || echo 0)
      la=\$(cut -d' ' -f1 /proc/loadavg)
      marker=\$( [ -f $W/state/stage2.done ] && echo DONE || echo running )
      echo \"\$(date +%H:%M) gff3=\$done_n bakta_procs=\$procs fails=\$fails load=\$la \$marker\"
      if [ \"\$marker\" = DONE ] || [ \"\$procs\" -eq 0 ]; then tail -6 $W/logs/stage2.out; echo STOPPOLL; fi
    " 2>&1)
  echo "$out"
  echo "$out" | grep -q STOPPOLL && exit 0
  sleep 300
done
echo "horizon reached"
