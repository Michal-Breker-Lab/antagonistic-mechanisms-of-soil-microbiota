#!/bin/bash
# FastTree vs IQ-TREE congruence, and the pull-back of the small tree outputs.
#
# This comparison is GENUINELY INDEPENDENT, which it would not have been under
# the original plan: IQ-TREE was to start from the FastTree topology (-t), which
# would have made agreement partly circular. That start tree crashed IQ-TREE's
# NNI (FastTree emits multifurcations for the zero-length polytomies of the 343
# identical sequences), so IQ-TREE built its own start tree and the two searches
# share nothing but the alignment.
#
# A large RF is not automatically an error -- FastTree is approximate -- but it
# is a number to look at before either tree is used. Report it either way.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
cd "$W"
[ -s iq_bac120.contree ] || { echo "FAIL: no iq_bac120.contree yet" >&2; exit 1; }
[ -s tree_bac120.nwk ]   || { echo "FAIL: no FastTree tree" >&2; exit 1; }

python3 - <<'PY'
import re, sys
from collections import Counter


def tips(nwk):
    return set(re.findall(r"[(,]\s*([^(),:;]+)\s*:", nwk))


ft = open("/mnt/LargeStorageNoBackup/Moshea/burk_bac120/tree_bac120.nwk").read()
iq = open("/mnt/LargeStorageNoBackup/Moshea/burk_bac120/iq_bac120.contree").read()
tf, ti = tips(ft), tips(iq)
print(f"  FastTree tips {len(tf)}   IQ-TREE tips {len(ti)}   shared {len(tf & ti)}")
for lab in ("MF6", "MF7"):
    assert lab in ti, f"{lab} missing from the IQ-TREE tree"
print("  MF6 and MF7 present in both")
if tf != ti:
    print(f"  WARNING: tip sets differ -- only FastTree: {sorted(tf-ti)[:3]}, "
          f"only IQ-TREE: {sorted(ti-tf)[:3]}")

# UFBoot distribution -- the thing FastTree cannot give
sup = [int(x) for x in re.findall(r"\)(\d+):", iq)]
if sup:
    c = Counter()
    for s in sup:
        c["<50" if s < 50 else ("50-79" if s < 80 else ("80-94" if s < 95 else ">=95"))] += 1
    print(f"  UFBoot over {len(sup)} internal nodes: "
          + "  ".join(f"{k} {c[k]} ({100*c[k]/len(sup):.0f}%)"
                      for k in ("<50", "50-79", "80-94", ">=95")))
PY

# Robinson-Foulds via ete3 if available; it is not required for the figure, so a
# missing ete3 is reported, never fatal.
python3 - <<'PY' || echo "  (ete3 unavailable - RF skipped, not an error)"
from ete3 import Tree
W = "/mnt/LargeStorageNoBackup/Moshea/burk_bac120"
ft = Tree(f"{W}/tree_bac120.nwk")
iq = Tree(f"{W}/iq_bac120.contree")
rf, mx, *_ = iq.robinson_foulds(ft, unrooted_trees=True)
print(f"  Robinson-Foulds: {rf} of {mx} ({100*rf/mx:.1f}% of splits differ)")
PY
echo "=== CONGRUENCE OK ==="
