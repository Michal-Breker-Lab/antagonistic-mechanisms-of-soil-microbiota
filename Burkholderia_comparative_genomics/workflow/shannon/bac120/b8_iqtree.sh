#!/bin/bash
# IQ-TREE full ML + 1000 UFBoot on the bac120 MSA -- THE PUBLISHED TREE.
#
# VERSION: 3.1.1, deliberately matched to the build that made this project's
# chromosome-1 and pC3 trees. Shannon also ships 3.0.1; using it would put a
# version difference between trees that get compared to each other.
#
# NO -t: the FastTree start tree was tried and CRASHES IQ-TREE --
#     "Assertion `node1->degree() == 3 && node2->degree() == 3' failed" in
#     getBestNNIForBran. FastTree emits MULTIFURCATING nodes for the zero-length
#     polytomies created by the 343 identical sequences, and IQ-TREE's NNI needs a
#     strictly bifurcating start tree. Arbitrarily resolving those polytomies would
#     inject a fabricated topology into the starting point, so IQ-TREE builds its
#     own start tree instead. The runtime argument for -t is gone anyway: this
#     alignment has 3,385 distinct site patterns against chr1's 677,041.
#     BONUS, and it is a real one: the FastTree-vs-IQ-TREE congruence check in
#     Task 6b Step 5 is now GENUINELY INDEPENDENT rather than partly circular.
# -mset: unrestricted -m MFP tests 200+ AA models, and ModelFinder ALONE was
#     measured at 10 h on this project's chr1 alignment. Restrict to the families
#     that actually win on bacterial marker data.
# --bcor 0.98 --nstep 10: the chr1 run lost ~24 h to UFBoot stalling at
#     correlation 0.983 against the 0.99 default, and .contree/.splits.nex are
#     written ONLY on completion -- cancelling loses every support value.
#     Loosen the target UP FRONT rather than discovering this again.
# --keep-ident IS REQUIRED, and the reason is not obvious. By default IQ-TREE
#     STRIPS sequences identical to another ("is ignored but added at the end")
#     and only then validates -t against the REDUCED alignment. This set is 54%
#     duplicates at 5,010 masked columns, so the FastTree start tree carries 259
#     taxa the reduced alignment no longer has, and the run aborts with
#     "Tree taxon X does not appear in the alignment". --keep-ident keeps every
#     sequence in the search, so the start tree matches AND all 791 genomes are
#     real tips rather than re-grafted after the fact.
#     Cost is negligible here: the alignment has only 3,385 distinct site
#     patterns (the chromosome-1 alignment had 677,041).
# NO -redo ANYWHERE: it discards the checkpoint. Re-running this script resumes.
set -euo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
IQ=/mnt/LargeStorageNoBackup/Moshea/envs/gtdbtk-2.7.2/bin/iqtree
cd "$W"; mkdir -p logs
CPUS=${CPUS:-64}

echo "=== preflight ==="
[ -s msa_bac120.faa ]  || { echo "FAIL: no MSA" >&2; exit 1; }
# tree_bac120.nwk is NOT a start tree any more, but b7 must still have run:
# it is the independent topology the congruence check compares against.
[ -s tree_bac120.nwk ] || { echo "FAIL: no FastTree tree -- run b7 first" >&2; exit 1; }
"$IQ" --version | head -1
for g in MF6 MF7; do
    grep -q "^>$g" msa_bac120.faa || { echo "FAIL: $g not in the MSA" >&2; exit 1; }
done
echo "  taxa: $(grep -c '^>' msa_bac120.faa)   threads: $CPUS"

"$IQ" -s msa_bac120.faa \
      -m MFP -mset LG,WAG,Q.pfam --mrate G4,R4 \
      -B 1000 --bcor 0.98 --nstep 10 --keep-ident \
      -T "$CPUS" --prefix iq_bac120

echo "=== verify ==="
for f in treefile contree splits.nex iqtree; do
    [ -s "iq_bac120.$f" ] || { echo "FAIL: iq_bac120.$f not written" >&2; exit 1; }
    printf "  iq_bac120.%-12s %s bytes\n" "$f" "$(stat -c%s "iq_bac120.$f")"
done
grep -E "Best-fit model|correlation coefficient|BEST SCORE" iq_bac120.log | tail -4 || true
python3 - <<'PY'
import re, sys
nwk = open('iq_bac120.contree').read()
tips = set(re.findall(r'[(,]\s*([^(),:]+):', nwk))
print(f"  tips: {len(tips)}")
if len(tips) != 791:
    print(f"  FAIL: expected 791 tips, got {len(tips)}"); sys.exit(1)
for g in ("MF6", "MF7"):
    if g not in tips:
        print(f"  FAIL: {g} absent"); sys.exit(1)
print("  MF6 and MF7 both present")
sup = [int(x) for x in re.findall(r'\)(\d+):', nwk)]
if sup:
    lo = sum(1 for s in sup if s < 50)
    print(f"  UFBoot nodes: {len(sup)}   <50: {lo} ({100*lo/len(sup):.1f}%)   "
          f">=95: {sum(1 for s in sup if s>=95)}")
PY
echo "=== IQTREE OK ==="
