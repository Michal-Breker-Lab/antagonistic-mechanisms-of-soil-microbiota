#!/bin/bash
# FastTree approximate ML on the bac120 MSA.
#
# ROLE: this is the SCAFFOLD, not the published tree. It gives IQ-TREE a start
# topology (most of the runtime saving) and an independent topology to check
# IQ-TREE against. Its supports are SH-like local supports, which are NOT
# bootstrap values and are NOT comparable to the chromosome-1 tree's UFBoot.
#
# DOUBLE PRECISION IS REQUIRED, not a preference: this set contains a near-clone
# block at ~0.0055 subs/site, and single-precision FastTree cannot resolve branch
# lengths at that scale -- it silently collapses them to zero. Shannon's build
# reports "Double precision" in its banner; assert it rather than trusting it.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
FT=/home/moshea/miniconda/envs/anvio-9/bin/FastTreeMP
cd "$W"; mkdir -p logs
export OMP_NUM_THREADS=${CPUS:-64}

echo "=== preflight ==="
[ -s msa_bac120.faa ] || { echo "FAIL: no MSA -- run b4 first" >&2; exit 1; }
BANNER=$("$FT" 2>&1 </dev/null | head -1 || true)
echo "  $BANNER"
case "$BANNER" in
    *"Double precision"*) ;;
    *) echo "FAIL: FastTree is not a double-precision build" >&2; exit 1;;
esac
TAXA=$(grep -c '^>' msa_bac120.faa)
echo "  taxa: $TAXA   threads: $OMP_NUM_THREADS"

echo "=== FastTree (LG + CAT) ==="
# -lg: the amino-acid model matching what ModelFinder picks for this data class.
# -gamma: rescale branch lengths under a gamma model so they are comparable to
#         IQ-TREE's, and report a real likelihood.
"$FT" -lg -gamma -log fasttree.log < msa_bac120.faa > tree_bac120.nwk

echo "=== verify ==="
[ -s tree_bac120.nwk ] || { echo "FAIL: no tree written" >&2; exit 1; }
python3 - <<'PY'
import re, sys
nwk = open('tree_bac120.nwk').read()
tips = re.findall(r'[(,]\s*([^(),:]+):', nwk)
print(f"  tips in tree: {len(tips)}")
missing = [g for g in ("MF6", "MF7") if g not in tips]
if missing:
    print(f"  FAIL: {missing} absent from the tree"); sys.exit(1)
print("  MF6 and MF7 both present")
# A tree of all-zero internal branches means precision was lost -- the exact
# failure double precision exists to prevent. Check rather than assume.
bl = [float(x) for x in re.findall(r':([0-9.eE-]+)', nwk)]
nz = sum(1 for b in bl if b > 0)
print(f"  branches: {len(bl)}   non-zero: {nz} ({100*nz/len(bl):.1f}%)")
print(f"  shortest non-zero: {min(b for b in bl if b > 0):.3g}")
PY
grep -E "^Total branch|ML_Lengths|Optimize all" fasttree.log | tail -3 || true
echo "=== FASTTREE OK ==="
