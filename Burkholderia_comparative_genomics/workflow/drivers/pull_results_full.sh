#!/usr/bin/env bash
# Pull the rebuild's small deliverables from Shannon to the Drive.
#
# Uses a single tar stream over ssh rather than rsync: rsync fails on every file
# over the DrvFs Drive mount (mkstemp/chgrp), and the original pull script hid
# that by sending stderr to /dev/null. A tar stream also avoids per-file
# round-trips, which matter on a 96%-full DrvFs mount.
#
# Genomes, annotations and Panaroo/mmseqs intermediates stay on Shannon: the
# Drive is a deliverables store with ~54 GB free, not a workspace.
set -uo pipefail
H=${H:-shannon}
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
D="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$D/tables" "$D/logs" "$D/trees"

pull() {   # pull <remote-dir> <local-dir> <find-expr...>
  local rdir="$1" ldir="$2"; shift 2
  echo "=== $rdir -> $(basename "$ldir") ==="
  ssh "$H" "cd '$rdir' 2>/dev/null && find . -maxdepth 1 -type f \\( $* \\) -print0 \
            | tar --null -cf - -T -" 2>/dev/null | tar xf - -C "$ldir" 2>/dev/null
  local rc=$?
  [ $rc -ne 0 ] && echo "  WARN: stream returned $rc"
  return 0
}

pull "$W/results" "$D/tables" "-name '*.tsv' -o -name '*.txt' -o -name '*.json'"
pull "$W/trees"   "$D/tables" "-name '*.treefile'"
pull "$W/trees"   "$D/logs"   "-name '*.iqtree'"
pull "$W/logs"    "$D/logs"   "-name '*.log'"

echo "=== pangenome matrices ==="
for m in chr1_strict_full c3_moderate_full; do
  if ssh "$H" "test -s $W/pangenome/$m/gene_presence_absence.Rtab"; then
    sz=$(ssh "$H" "stat -c%s $W/pangenome/$m/gene_presence_absence.Rtab")
    # chr1 was 84 MB at 306 genomes; at 771 it will be far larger. Drive has
    # ~54 GB free, so anything over 500 MB stays on Shannon by default.
    if [ "$sz" -gt 524288000 ] && [ "${FORCE_RTAB:-0}" != "1" ]; then
      echo "  $m Rtab is $((sz/1048576)) MB - left on Shannon (FORCE_RTAB=1 to override)"
    else
      ssh "$H" "cat $W/pangenome/$m/gene_presence_absence.Rtab" \
        > "$D/tables/${m}_gene_presence_absence.Rtab"
      echo "  ${m}_gene_presence_absence.Rtab ($((sz/1048576)) MB)"
    fi
    ssh "$H" "cat $W/pangenome/$m/summary_statistics.txt" \
      > "$D/tables/${m}_summary_statistics.txt" 2>/dev/null
  else
    echo "  $m MISSING"
  fi
done

echo
echo "=== what landed ==="
du -sh "$D/tables" "$D/logs"
echo "Drive free: $(df -h "$D" | tail -1 | awk '{print $4}')"
