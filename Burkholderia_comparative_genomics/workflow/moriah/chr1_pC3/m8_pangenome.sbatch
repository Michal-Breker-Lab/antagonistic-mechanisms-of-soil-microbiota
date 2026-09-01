#!/bin/bash
#SBATCH --job-name=c3_pang
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/m8_pang.%j.out
#SBATCH --error=logs/m8_pang.%j.err
# NO --partition. --mem and --time are supplied per-arm at submit time.
#
# Parameterised pangenome arm. Required env (via --export):
#   SET   = pc3 | chr1
#   TOOL  = panaroo | ppanggolin
#   MODE  = panaroo clean-mode (strict|moderate|sensitive)  OR  ppanggolin identity (0.80|0.60)
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
PY="$M/bin/python"
THREADS=${SLURM_CPUS_PER_TASK:-32}
: "${SET:?SET not set}"; : "${TOOL:?TOOL not set}"; : "${MODE:?MODE not set}"
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/pang.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

TAG="${SET}_${TOOL}_${MODE//./}"
GFFDIR="$W/gff_${SET}"
OUT="$W/pangenome_${TAG}"

echo "=== preflight (RUN the tools, and CHECK the inputs) ==="
LAB="$W/results/pangenome_${SET}_labels.tsv"
[ -s "$LAB" ] || { echo "PREFLIGHT FAIL: missing $LAB" >&2; exit 1; }
echo "  input ok: $LAB ($(( $(wc -l < "$LAB") - 1 )) replicons)"
source "$M/etc/profile.d/conda.sh"
if [ "$TOOL" = "panaroo" ]; then
    conda activate panaroo && panaroo --version
else
    # ACTIVATE the env -- an absolute path to the binary is not enough.
    # PPanGGOLiN shells out to `mmseqs` by bare name, so its env bin must be on
    # PATH or clustering dies with "Command 'mmseqs' not found" a few minutes in.
    conda activate ppanggolin
    ppanggolin --version
    command -v mmseqs >/dev/null || {
        echo "PREFLIGHT FAIL: mmseqs not on PATH after activating ppanggolin" >&2
        exit 1; }
    echo "  mmseqs on PATH: $(command -v mmseqs) ($(mmseqs version))"
fi
echo "SET=$SET TOOL=$TOOL MODE=$MODE THREADS=$THREADS TMPDIR=$TMPDIR"

# --- subset GFF3s once per SET (shared by all arms of that set) -------------
EXPECT=$(( $(wc -l < "$LAB") - 1 ))
if [ ! -d "$GFFDIR" ] || [ "$(ls "$GFFDIR"/*.gff3 2>/dev/null | wc -l)" -ne "$EXPECT" ]; then
    echo "=== subsetting GFF3s for $SET ==="
    "$PY" "$W/s27_subset_gff.py" --labels "$LAB" --annot "$W/annot" --outdir "$GFFDIR"
else
    echo "=== GFF3s for $SET already present ($EXPECT) ==="
fi
N=$(ls "$GFFDIR"/*.gff3 | wc -l)
[ "$N" -eq "$EXPECT" ] || { echo "FAIL: $N GFF3s, expected $EXPECT" >&2; exit 1; }

# --- run the arm ------------------------------------------------------------
if [ "$TOOL" = "panaroo" ]; then
    echo "=== Panaroo ($MODE) on $N genomes ==="
    panaroo -i "$GFFDIR"/*.gff3 -o "$OUT" --clean-mode "$MODE" \
            -t "$THREADS" --remove-invalid-genes --aligner mafft -a core
else
    echo "=== PPanGGOLiN (identity $MODE) on $N genomes ==="
    ls "$GFFDIR"/*.gff3 | while read -r f; do
        printf "%s\t%s\n" "$(basename "$f" .gff3)" "$f"
    done > "$TMPDIR/ppg_input.tsv"
    ppanggolin workflow --anno "$TMPDIR/ppg_input.tsv" \
        --output "$OUT" --identity "$MODE" --coverage 0.8 --cpu "$THREADS" --force
fi

echo "=== DONE $TAG ==="
exit 0
