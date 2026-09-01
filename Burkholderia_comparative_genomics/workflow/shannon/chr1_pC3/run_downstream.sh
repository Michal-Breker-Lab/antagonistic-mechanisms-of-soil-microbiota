#!/usr/bin/env bash
# Chain every downstream stage once Bakta annotation finishes, so the pipeline
# keeps running unattended. Each stage logs separately and the chain stops at the
# first hard failure rather than feeding empty inputs forward.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
cd $W
L=$W/logs/downstream.log
say () { echo "[$(date '+%F %T')] $*" | tee -a $L; }

say "waiting for annotation to finish..."
while ! grep -q ANNOT_DONE $W/logs/s6_annotate.log 2>/dev/null; do sleep 120; done
n=$(ls $W/annot/*/*.gff3 2>/dev/null | wc -l)
f=$(grep -c '^FAIL' $W/logs/s6_annotate.log 2>/dev/null || echo 0)
say "annotation complete: $n genomes, $f failures"

say "STAGE 4a - replicon typing"
python3 -u s4a_typing.py > logs/s4a.log 2>&1 || { say "4a FAILED"; exit 1; }
say "  $(grep -c . logs/s4a.log) lines; $(grep -m1 'secondary_large' logs/s4a.log)"

say "STAGE 4b - ortholog clustering"
bash s4b_coherence.sh > logs/s4b.log 2>&1 || { say "4b FAILED"; exit 1; }
say "  $(grep -m1 'clusters:' logs/s4b.log)"

say "STAGE 4c - c3 coherence + calls"
python3 -u s4c_identity.py > logs/s4c.log 2>&1 || { say "4c FAILED"; exit 1; }
say "  $(grep -m1 'replicons called c3' logs/s4c.log)"

say "STAGE 7+8 - pangenomes (this is the long one)"
bash s78_pangenomes.sh > logs/s78.log 2>&1
say "  panaroo done; chr1 core: $(awk -F'\t' '/^Core genes/{print $3}' \
     $W/pangenome/chr1_strict/summary_statistics.txt 2>/dev/null)"
say "  c3 core: $(awk -F'\t' '/^Core genes/{print $3}' \
     $W/pangenome/c3_moderate/summary_statistics.txt 2>/dev/null)"

say "STAGE 8+9 - trees"
bash s89_trees.sh > logs/s89.log 2>&1
say "  trees: $(ls $W/trees/*.treefile 2>/dev/null | wc -l)"

say "STAGE 11 - phylogenetic statistics"
if [ -s $W/trees/chr1_core.treefile ]; then
  W=$W ~/miniconda/envs/R_analysis/bin/Rscript s11_phylo_stats.R > logs/s11.log 2>&1 \
    && say "  stats written" || say "  stats FAILED (see logs/s11.log)"
else
  say "  no species tree; skipping stats"
fi

say "DOWNSTREAM_ALL_DONE"
