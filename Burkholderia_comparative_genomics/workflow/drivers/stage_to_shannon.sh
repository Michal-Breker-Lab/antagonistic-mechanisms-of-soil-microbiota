#!/usr/bin/env bash
# Copy the rebuild scripts and the MF7 assembly to Shannon.
#
# scp, not rsync: rsync dies on every file over the DrvFs Drive mount
# (mkstemp/chgrp failures), and piping it to tail hides the failure.
set -euo pipefail
H=${H:-shannon}
W=${W:-/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3}
HERE="$(cd "$(dirname "$0")" && pwd)"
DRIVE="$(cd "$HERE/.." && pwd)"
MF7="/mnt/g/My Drive/Moshe/Collaborations/Tzila_Results/MF7/IMG_2546825540/IMG Data/MF7.fna"

echo "=== staging scripts to $H:$W ==="
for f in s28_clone_labels.py s29_null_stratified.py s30_mf7_replicon_assign.py \
         clone_collapse.py s6b_annotate_rest.sh run_full_rebuild.sh watch_shannon.sh; do
  [ -f "$HERE/$f" ] || { echo "missing $f" >&2; exit 1; }
  scp -q "$HERE/$f" "$H:$W/$f"
  echo "  $f"
done

echo "=== staging MF7 ==="
if [ ! -f "$MF7" ]; then echo "MF7.fna not found at: $MF7" >&2; exit 1; fi
scp -q "$MF7" "$H:$W/MF7.fna"
ssh "$H" "ls -l $W/MF7.fna; grep -c '>' $W/MF7.fna"

echo "=== verifying ==="
ssh "$H" "cd $W && python3 -c 'import ast;[ast.parse(open(f).read()) for f in [\"s28_clone_labels.py\",\"s29_null_stratified.py\",\"s30_mf7_replicon_assign.py\",\"clone_collapse.py\"]];print(\"python OK\")' && bash -n run_full_rebuild.sh && bash -n s6b_annotate_rest.sh && echo 'bash OK'"

echo
echo "Staged. Self-test the null machinery on Shannon before relying on it:"
echo "  ssh $H 'cd $W && python3 s29_null_stratified.py --selftest'"
echo "Then, when check_shannon.sh says FREE:"
echo "  ssh $H 'cd $W && nohup bash run_full_rebuild.sh --stage 0 > logs/stage0.out 2>&1'"
