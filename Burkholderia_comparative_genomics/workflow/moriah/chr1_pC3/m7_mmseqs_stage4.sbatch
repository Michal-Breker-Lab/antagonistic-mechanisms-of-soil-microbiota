#!/bin/bash
#SBATCH --job-name=c3_stage4
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=06:00:00
#SBATCH --output=logs/m7_stage4.%j.out
#SBATCH --error=logs/m7_stage4.%j.err
# NO --partition: the -k variants are account-restricted.
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
PY="$M/bin/python"
MMSEQS="$M/envs/hgt_recon/bin/mmseqs"
THREADS=${SLURM_CPUS_PER_TASK:-32}
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/stage4.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight (RUN the tools, and CHECK the inputs) ==="
"$MMSEQS" version
"$PY" -c "import sys; print('python', sys.version.split()[0])"
for f in "$W/results/replicon_types.tsv" "$W/results/pc3_training_labels.tsv"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing input $f" >&2; exit 1; }
    echo "  input ok: $f ($(wc -l < "$f") lines)"
done
N=$(ls -d "$W"/annot/*/ 2>/dev/null | wc -l)
[ "$N" -ge 773 ] || { echo "PREFLIGHT FAIL: only $N annot dirs" >&2; exit 1; }
echo "  input ok: $N annotation directories"
echo "TMPDIR=$TMPDIR  THREADS=$THREADS"

# --- Stage 4a: collect proteins from every large secondary replicon ---------
echo "=== Stage 4a: secondary-replicon proteins ==="
"$PY" "$W/s31_secondary_proteins.py" --types "$W/results/replicon_types.tsv" \
    --annot "$W/annot" --out "$W/results/secondary_proteins.faa"

# --- Stage 4b: MMseqs2 easy-cluster, 50% identity / 80% coverage ------------
# Parameters taken verbatim from report_methods_draft.md.
echo "=== Stage 4b: mmseqs easy-cluster (50% id, 80% cov) ==="
"$MMSEQS" easy-cluster "$W/results/secondary_proteins.faa" \
    "$W/results/secondary_clu" "$TMPDIR/mm" \
    --min-seq-id 0.5 -c 0.8 --threads "$THREADS"
wc -l "$W/results/secondary_clu_cluster.tsv"

# --- Stage 4c: residence, diagnostics, classification ----------------------
echo "=== Stage 4c: diagnostics (frac>=0.90, n>=10) ==="
"$PY" "$W/s26_pc3_diagnostics.py" --clusters "$W/results/secondary_clu_cluster.tsv" \
    --labels "$W/results/pc3_training_labels.tsv" \
    --types "$W/results/replicon_types.tsv" \
    --frac 0.90 --min-genomes 10 --outdir "$W/results"

echo "=== D9c sensitivity at 0.95 ==="
"$PY" "$W/s26_pc3_diagnostics.py" --clusters "$W/results/secondary_clu_cluster.tsv" \
    --labels "$W/results/pc3_training_labels.tsv" \
    --types "$W/results/replicon_types.tsv" \
    --frac 0.95 --min-genomes 10 --outdir "$W/results" --tag frac95

echo "=== DONE ==="
exit 0
