#!/usr/bin/env bash
# Retry Shannon until sshd answers, then report stage-1 state and (re)launch the
# chain if it is not running. Shannon's sshd refuses connections outright when
# the box is in a runaway state, so a single failed ssh means nothing.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
for i in $(seq 1 240); do            # up to ~60 x 60 s
  if out=$(timeout 45 ssh -o BatchMode=yes -o ConnectTimeout=20 shannon "
        head -1 /proc/loadavg
        echo '--- stage1 ---'; tail -3 $W/logs/stage1.out 2>/dev/null
        echo '--- markers ---'; ls $W/state/ 2>/dev/null | tr '\n' ' '
        echo; echo '--- chain ---'
        if ps -eo cmd | grep -qE '[c]hain_after_stage1'; then echo running; else
          cd $W && setsid nohup bash chain_after_stage1.sh >/dev/null 2>&1 </dev/null &
          sleep 2; echo '(re)launched'
        fi
        tail -3 $W/logs/chain1.log 2>/dev/null
     " 2>&1); then
    echo "=== reached Shannon on attempt $i ($(date +%T)) ==="
    echo "$out"
    exit 0
  fi
  sleep 60
done
echo "never reachable in 60 attempts"
exit 1
