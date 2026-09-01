#!/bin/bash
#SBATCH --job-name=c3_iqtree
#SBATCH --partition=glacier
# --partition is normally omitted (house rule: the -k variants are
# account-restricted and SLURM picks by account). It is set explicitly here
# because the site job_submit plugin stamps Features=cpuonly but sometimes
# leaves the full partition list (glacier,catfish,salmon,goldfish,dogfish); the
# GPU partitions cannot satisfy cpuonly, so SLURM 20.11 rejects the whole
# request as BadConstraints instead of falling back. glacier is the non-`-k`
# partition SLURM selected for every job that ran today.
#SBATCH --output=logs/m10_iqtree.%j.out
#SBATCH --error=logs/m10_iqtree.%j.err
set -euo pipefail
# Required env: SET = pc3 | chr1 ; optional START = newick starting tree,
#               MODEL (default GTR+F+I+R5)
W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
IQ="$M/envs/anvio-9/bin/iqtree3"
T=${SLURM_CPUS_PER_TASK:-32}
: "${SET:?SET not set}"
cd "$W"
export TMPDIR=$(mktemp -d "$W/tmp/iq.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT

# ALN may be overridden to use a taxon-reduced alignment (see D13).
if [ -n "${ALN:-}" ]; then
    :
else
  case "$SET" in
    pc3)  ALN="$W/pangenome_pc3_panaroo_moderate/core_gene_alignment.aln" ;;
    chr1) ALN="$W/pangenome_chr1_panaroo_strict/core_gene_alignment.aln" ;;
    *) echo "bad SET=$SET" >&2; exit 1 ;;
  esac
fi

echo "=== preflight ==="
[ -s "$ALN" ] || { echo "PREFLIGHT FAIL: missing $ALN" >&2; exit 1; }
echo "  input ok: $ALN ($(grep -c '^>' "$ALN") taxa)"
"$IQ" --version | head -1
EXTRA=""
if [ -n "${START:-}" ]; then
    [ -s "$START" ] || { echo "PREFLIGHT FAIL: missing start tree $START" >&2; exit 1; }
    echo "  start tree: $START"
    EXTRA="-t $START"
fi

# Substitution model is NOT re-searched with ModelFinder. It is taken from the
# ORIGINAL run's own ModelFinder output, retained on the Drive in
# logs/c3_core.iqtree and logs/chr1_core.iqtree:
#
#   pC3  core (140 seq x 52,834 sites)  -> GTR+F+I+R5
#   chr1 core (304 seq x 490,720 sites) -> GTR+F+R10
#
# Both selected GTR+F. Rate heterogeneity is strongly multimodal here -- chr1's
# fitted R10 puts 46% of sites at rate 0.012 with a thin tail out to 13.8 -- so
# a 4-category model is badly under-parameterised. GTR+F+I+R5 is used for both
# trees: it is the original's exact selection for pC3, and for chr1 it captures
# the large near-invariant class with an explicit +I rather than spending
# categories on it, at roughly half the per-site cost of R10.
#
# This reuses the original ModelFinder result rather than repeating the search,
# so it is faithful AND free. State in Methods that the model was inherited, not
# re-selected.
MODEL=${MODEL:-GTR+F+I+R5}
echo "=== IQ-TREE ($SET), $MODEL, 1000 UFBoot ==="
"$IQ" -s "$ALN" -m "$MODEL" -B 1000 -T "$T" $EXTRA \
      --prefix "$W/results/tree_${SET}" -redo

echo "=== DONE ==="
exit 0
