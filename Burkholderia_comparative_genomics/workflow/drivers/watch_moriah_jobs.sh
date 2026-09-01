#!/bin/bash
# Watch ALL of this project's Moriah jobs until none remain, then report each
# one's terminal state. Requires an SSH_OK sentinel before concluding a job left
# the queue: a dropped VPN otherwise looks identical to completion, which is how
# an earlier watcher reported a false "LEFT QUEUE".
set -uo pipefail
JOBS="$*"
[ -n "$JOBS" ] || { echo "usage: $0 <jobid> [jobid...]" >&2; exit 1; }
W=/sci/backup/ofinkel/moshea/burkholderia_c3

while :; do
    OUT=$(ssh -o ConnectTimeout=20 -o BatchMode=yes moriah \
          "bash -lc 'echo SSH_OK; squeue -u moshea -h -o \"%i %T\"'" 2>/dev/null)
    if ! grep -q '^SSH_OK$' <<<"$OUT"; then
        echo "$(date +%H:%M) unreachable (VPN/ssh) - retrying"
        sleep 180; continue
    fi
    RUNNING=""
    for j in $JOBS; do
        grep -q "^$j " <<<"$OUT" && RUNNING="$RUNNING $j"
    done
    if [ -z "$RUNNING" ]; then
        echo "$(date +%H:%M) all watched jobs have left the queue"
        break
    fi
    echo "$(date +%H:%M) still running:$RUNNING"
    sleep 300
done

ssh -o BatchMode=yes moriah "bash -lc '
for j in $JOBS; do
    echo \"--- \$j ---\"
    sacct -j \$j --format=JobID%14,JobName%14,State,Elapsed,ExitCode -n | head -2
done'" 2>&1
