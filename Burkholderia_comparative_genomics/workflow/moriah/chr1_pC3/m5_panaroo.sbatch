#!/bin/bash
#SBATCH --job-name=c3_panaroo
#SBATCH --cpus-per-task=32
#SBATCH --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=logs/m5_panaroo.%j.out
#SBATCH --error=logs/m5_panaroo.%j.err
# NO --partition: the -k variants are account-restricted; SLURM picks by account.
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
PY="$M/bin/python"
THREADS=${SLURM_CPUS_PER_TASK:-32}
MODE=${PANAROO_MODE:-moderate}
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/panaroo.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

source "$M/etc/profile.d/conda.sh"
conda activate panaroo

echo "=== preflight (RUN the tools, and CHECK the inputs) ==="
panaroo --version
"$PY" -c "import sys; print('python', sys.version.split()[0])"
for f in "$W/results/pc3_training_labels.tsv" "$W/results/replicon_types.tsv"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing input $f" >&2; exit 1; }
    echo "  input ok: $f ($(wc -l < "$f") lines)"
done
N=$(ls "$W"/training_gff/*.gff3 2>/dev/null | wc -l)
[ "$N" -eq 136 ] || { echo "PREFLIGHT FAIL: expected 136 training GFF3s, found $N" >&2; exit 1; }
echo "  input ok: 136 training GFF3s"
echo "TMPDIR=$TMPDIR  THREADS=$THREADS  MODE=$MODE"

OUT="$W/panaroo_training_${MODE}"
echo "=== Panaroo (${MODE}) on 136 secondary-replicon subsets ==="
panaroo -i "$W"/training_gff/*.gff3 -o "$OUT" --clean-mode "$MODE" \
        -t "$THREADS" --remove-invalid-genes

echo "=== derive pC3 diagnostics (D9, frac>=0.90, n>=10) ==="
"$PY" "$W/s26_pc3_diagnostics.py" \
    --presence "$OUT/gene_presence_absence.csv" \
    --gene-data "$OUT/gene_data.csv" \
    --labels "$W/results/pc3_training_labels.tsv" \
    --types "$W/results/replicon_types.tsv" \
    --frac 0.90 --min-genomes 10 --outdir "$W/results"

echo "=== D9c sensitivity: same derivation at 0.95 ==="
"$PY" "$W/s26_pc3_diagnostics.py" \
    --presence "$OUT/gene_presence_absence.csv" \
    --gene-data "$OUT/gene_data.csv" \
    --labels "$W/results/pc3_training_labels.tsv" \
    --types "$W/results/replicon_types.tsv" \
    --frac 0.95 --min-genomes 10 --outdir "$W/results" --tag frac95

echo "=== DONE ==="
exit 0
