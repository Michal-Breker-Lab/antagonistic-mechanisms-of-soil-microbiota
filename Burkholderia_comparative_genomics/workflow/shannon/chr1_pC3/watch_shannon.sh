#!/usr/bin/env bash
# Detached watcher: writes SHANNON_FREE.flag once the machine is genuinely idle.
# Runs ON Shannon. Does NOT launch anything -- the launch stays a human decision.
set -uo pipefail
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
NEED=${NEED:-60}          # free cores required
STREAK=${STREAK:-3}       # consecutive samples before declaring free
INTERVAL=${INTERVAL:-300}
FLAG="$W/logs/SHANNON_FREE.flag"
LOG="$W/logs/watch_shannon.log"
mkdir -p "$W/logs"
echo $$ > "$W/logs/watch_shannon.pid"   # stop by PIDFILE, never pkill -f

cores=$(nproc); hits=0
while true; do
  read -r _ u1 n1 s1 i1 w1 _ < /proc/stat; sleep 5
  read -r _ u2 n2 s2 i2 w2 _ < /proc/stat
  dt=$(( (u2-u1)+(n2-n1)+(s2-s1)+(i2-i1)+(w2-w1) ))
  di=$(( (i2-i1)+(w2-w1) ))
  free=$(awk -v di="$di" -v dt="$dt" -v c="$cores" 'BEGIN{printf "%.0f",(dt>0?di/dt*c:0)}')
  if [ "$free" -ge "$NEED" ]; then hits=$((hits+1)); else hits=0; fi
  echo "$(date +%F' '%T) free=${free}/${cores} streak=${hits}/${STREAK}" >> "$LOG"
  if [ "$hits" -ge "$STREAK" ]; then
    date +%F' '%T > "$FLAG"
    echo "$(date +%F' '%T) FREE - flag written, watcher exiting" >> "$LOG"
    exit 0
  fi
  sleep "$INTERVAL"
done
