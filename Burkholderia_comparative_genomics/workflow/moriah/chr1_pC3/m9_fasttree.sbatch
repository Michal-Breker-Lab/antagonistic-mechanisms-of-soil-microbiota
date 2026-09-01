#!/bin/bash
#SBATCH --job-name=c3_fasttree
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=glacier
# --partition is normally omitted (house rule: the -k variants are
# account-restricted and SLURM picks by account). It is set explicitly here
# because the site job_submit plugin stamps Features=cpuonly but sometimes
# leaves the full partition list (glacier,catfish,salmon,goldfish,dogfish); the
# GPU partitions cannot satisfy cpuonly, so SLURM 20.11 rejects the whole
# request as BadConstraints instead of falling back. glacier is the non-`-k`
# partition SLURM selected for every job that ran today.
#SBATCH --output=logs/m9_fasttree.%j.out
#SBATCH --error=logs/m9_fasttree.%j.err
set -euo pipefail
W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
cd "$W"
export TMPDIR=$(mktemp -d "$W/tmp/ft.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

ALN="$W/pangenome_chr1_panaroo_strict/core_gene_alignment.aln"
echo "=== preflight ==="
[ -s "$ALN" ] || { echo "PREFLIGHT FAIL: missing $ALN" >&2; exit 1; }
echo "  input ok: $ALN ($(grep -c '^>' "$ALN") taxa)"
"$M/envs/anvio-9/bin/FastTree" -help 2>&1 | head -1

echo "=== FastTree (GTR+CAT) on chr1 core alignment ==="
# Starting topology for the IQ-TREE search, and an independent check on it.
"$M/envs/anvio-9/bin/FastTree" -nt -gtr -nosupport \
    "$ALN" > "$W/results/chr1_fasttree.nwk"
echo "tips: $(tr ',' '\n' < "$W/results/chr1_fasttree.nwk" | grep -c ':')"
echo "=== DONE ==="
exit 0
