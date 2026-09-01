#!/usr/bin/env bash
# Retry Shannon with a generous connect timeout. The box answers ICMP but its
# sshd cannot get scheduled to send a banner when loaded, so short timeouts
# guarantee failure regardless of whether it is recoverable. 120 s each.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
for i in $(seq 1 40); do
  if out=$(timeout 150 ssh -o BatchMode=yes -o ConnectTimeout=120 \
             -o ServerAliveInterval=20 shannon "
        date +%H:%M:%S
        echo \"load : \$(cut -d' ' -f1-3 /proc/loadavg)\"
        echo \"gff3 : \$(ls $W/annot/*/*.gff3 2>/dev/null | wc -l) / 773\"
        echo \"bakta: \$(ps -eo cmd | grep -cE '[b]akta')\"
        echo \"mine : \$(ps -u \$(id -u) -o cmd | grep -cE '[b]akta')\"
        echo \"other: \$(ps -eo user,cmd --no-headers | grep -E '[b]akta|python3|diamond' | grep -vc \"^\$(id -un)\")\"
        echo \"top  :\"; ps -eo user,pcpu,comm --sort=-pcpu --no-headers | head -5
        echo \"stage2: \$( [ -f $W/state/stage2.done ] && echo DONE || echo running )\"
     " 2>&1); then
    echo "=== reached on attempt $i ($(date +%T)) ==="; echo "$out"; exit 0
  fi
  echo "attempt $i failed ($(date +%T))"
  sleep 90
done
echo "unreachable after 40 attempts"
