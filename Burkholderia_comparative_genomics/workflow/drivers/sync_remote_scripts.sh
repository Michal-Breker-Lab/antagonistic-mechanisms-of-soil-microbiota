#!/usr/bin/env bash
# Mirror a remote host's project scripts onto the Drive.
#
# Required by the "Remote scripts MUST be mirrored to the Drive" rule in
# CLAUDE.md, added 2026-08-20 after Shannon went down and stranded an entire
# pipeline that amounted to a few hundred KB of text.
#
# Scripts and small irreplaceable inputs only. Genomes, annotations, databases
# and intermediates stay remote: the Drive is a deliverables store, and those
# regenerate.
#
# Transport notes learned the hard way:
#   - rsync fails on every file over the DrvFs Drive mount (mkstemp/chgrp).
#   - Moriah's sshd has NO sftp subsystem, so modern scp fails with
#     "subsystem request failed"; scp -O (legacy protocol) works.
#   - A tar stream over ssh sidesteps both and is one round-trip.
set -uo pipefail

HOST=${1:-}
case "$HOST" in
  shannon) REMOTE=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3 ;;
  moriah)  REMOTE=/sci/backup/ofinkel/moshea/burkholderia_c3 ;;
  "")      echo "usage: $0 <shannon|moriah> [remote_dir]" >&2; exit 2 ;;
  *)       REMOTE=${2:-}; [ -n "$REMOTE" ] || { echo "give remote dir for $HOST" >&2; exit 2; } ;;
esac
[ -n "${2:-}" ] && REMOTE=$2

D="$(cd "$(dirname "$0")" && pwd)/$HOST"
mkdir -p "$D"

echo "=== $HOST:$REMOTE -> $D ==="
if ! timeout 40 ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" true 2>/dev/null; then
  echo "  $HOST unreachable - nothing mirrored" >&2
  exit 1
fi

# Top-level scripts and small inputs. -maxdepth 1 keeps genomes/ and annot/ out.
timeout 600 ssh -o BatchMode=yes "$HOST" "
  cd '$REMOTE' 2>/dev/null || exit 1
  find . -maxdepth 1 -type f \\( \
      -name '*.sh' -o -name '*.sbatch' -o -name '*.py' -o -name '*.R' \
      -o -name '*.faa' -o -name '*.smk' -o -name 'Snakefile' \
      -o -name '*.yaml' -o -name '*.yml' \\) -size -5M -print0 \
  | tar --null -cf - -T - 2>/dev/null" | tar xf - -C "$D" 2>/dev/null

# sbatch/ subdir if the host uses one (Moriah does)
timeout 300 ssh -o BatchMode=yes "$HOST" "
  cd '$REMOTE/sbatch' 2>/dev/null && tar cf - . 2>/dev/null" \
  | tar xf - -C "$D" 2>/dev/null

n=$(find "$D" -type f ! -name 'desktop.ini' | wc -l)
echo "  mirrored $n files"
find "$D" -type f ! -name 'desktop.ini' -printf '    %f\n' 2>/dev/null | sort | head -40
du -sh "$D"
