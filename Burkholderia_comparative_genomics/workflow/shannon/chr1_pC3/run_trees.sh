#!/usr/bin/env bash
# Trees, run in parallel (48 threads each of 128 cores).
# Model search restricted to the GTR family (-mset GTR): on a 490 kb bacterial
# core-genome alignment ModelFinder's unrestricted search costs hours and
# effectively always lands on GTR+F+R/I+G anyway. Rate heterogeneity is still
# selected freely.
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
IQ=$W/envs/burk/bin/iqtree3
mkdir -p $W/trees
$IQ -s $W/pangenome/chr1_strict/core_gene_alignment.aln --prefix $W/trees/chr1_core \
    -m MFP -mset GTR -B 1000 -T 48 --seqtype DNA -redo > $W/logs/iqtree_chr1.log 2>&1 &
P1=$!
$IQ -s $W/pangenome/c3_moderate/core_gene_alignment.aln --prefix $W/trees/c3_core \
    -m MFP -mset GTR -B 1000 -T 40 --seqtype DNA -redo > $W/logs/iqtree_c3.log 2>&1 &
P2=$!
wait $P1; echo "chr1 exit=$?"
wait $P2; echo "c3 exit=$?"
for p in chr1_core c3_core; do
  [ -s $W/trees/$p.treefile ] && echo "$p OK: $(grep -m1 'Best-fit model' $W/trees/$p.iqtree)"
done
echo TREES_DONE
