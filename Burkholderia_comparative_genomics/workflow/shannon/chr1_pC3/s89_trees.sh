#!/usr/bin/env bash
# Stage 8 + 9 - species tree (chromosome-1 core) and c3 tree (c3 core).
#
# The species tree is built from chromosome 1 ONLY so that c3 presence cannot
# influence the topology onto which c3 presence is later mapped.
#
# Model selection is left to ModelFinder rather than hard-coded, so the
# substitution model is chosen by BIC from the data instead of by my guess.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
IQ=$W/envs/burk/bin/iqtree3
mkdir -p $W/trees $W/logs

build () {   # $1=alignment  $2=prefix  $3=label
  local aln=$1 pre=$2 lab=$3
  if [ ! -s "$aln" ]; then echo "SKIP $lab: no alignment at $aln"; return 1; fi
  local n=$(grep -c '^>' "$aln")
  local L=$(python3 -c "
import sys
seq=[];cur=[]
for l in open('$aln'):
    if l.startswith('>'):
        if cur: seq.append(''.join(cur)); cur=[]
    else: cur.append(l.strip())
if cur: seq.append(''.join(cur))
print(len(seq[0]) if seq else 0)")
  echo "=== $lab: $n taxa, $L bp ==="
  date
  $IQ -s "$aln" --prefix "$pre" -m MFP -B 1000 -T 48 --seqtype DNA -redo \
      > $W/logs/iqtree_$(basename $pre).log 2>&1
  echo "exit=$?"
  [ -s "$pre.treefile" ] && echo "tree written: $pre.treefile" \
      && grep -m1 'Best-fit model' "$pre.iqtree" 2>/dev/null
  date
}

build $W/pangenome/chr1_strict/core_gene_alignment.aln  $W/trees/chr1_core "species tree (chr1 core)"
build $W/pangenome/c3_moderate/core_gene_alignment.aln  $W/trees/c3_core   "c3 tree (c3 core)"

# ---- report c3 core size: decides whether the c3 tree is meaningful at all ----
echo "=== c3 pangenome core size check ==="
if [ -s $W/pangenome/c3_moderate/summary_statistics.txt ]; then
  cat $W/pangenome/c3_moderate/summary_statistics.txt
  core=$(awk -F'\t' '/^Core genes/{print $3}' $W/pangenome/c3_moderate/summary_statistics.txt)
  echo "c3 core genes: ${core:-unknown}"
  if [ -n "${core:-}" ] && [ "${core:-0}" -lt 50 ] 2>/dev/null; then
    echo "WARNING: c3 core < 50 genes. Per the pre-declared fallback (plan D7),"
    echo "the c3 phylogeny should be replaced by a gene-content dendrogram and"
    echo "labelled as such, NOT presented as a phylogeny."
  fi
fi
echo "TREES_DONE"
