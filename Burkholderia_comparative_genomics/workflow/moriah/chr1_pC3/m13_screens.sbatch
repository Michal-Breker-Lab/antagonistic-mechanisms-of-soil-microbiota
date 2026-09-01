#!/bin/bash
#SBATCH --job-name=c3_screen
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0-10:00:00
#SBATCH --output=logs/m13_screen.%j.out
#SBATCH --error=logs/m13_screen.%j.err
# NO --partition (the -k variants are account-restricted).
#
# The three fig7/fig8 screens over Set B (D15: 6 genomes incl. MF7, ANI >= 96%).
# Flags reproduce TOOLS.md; only -cpu/--cpus differ where the core count differs.
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
IPS=/sci/backup/ofinkel/moshea/Tools/interproscan-5.76-107.0
IN="$W/setb_pc3"; OUT="$W/screens"; DB="$W/databases/antismash_db"
MODELS="$W/macsy_models"
T=${SLURM_CPUS_PER_TASK:-32}
cd "$W"; mkdir -p "$OUT"

# TMPDIR on the lab FS: compute-node /tmp is per-user quota'd and InterProScan
# writes large temporaries (typed memory: moriah-tmp-quota-on-array-jobs).
export TMPDIR=$(mktemp -d "$W/tmp/screen.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight: CHECK THE INPUTS ==="
MEM=$(tail -n +2 "$W/setB_members.tsv" | cut -f1)
N=$(wc -l <<<"$MEM")
[ "$N" -eq 6 ] || { echo "PREFLIGHT FAIL: $N Set B members, expected 6" >&2; exit 1; }
for g in $MEM; do
    for e in fna gbk faa; do
        [ -s "$IN/${g}_pC3.$e" ] || { echo "PREFLIGHT FAIL: missing $IN/${g}_pC3.$e" >&2; exit 1; }
    done
done
[ -s "$IN/setB_pC3_all.faa" ] || { echo "PREFLIGHT FAIL: missing setB_pC3_all.faa" >&2; exit 1; }
[ -x "$IPS/interproscan.sh" ] || { echo "PREFLIGHT FAIL: no interproscan.sh" >&2; exit 1; }
[ -d "$DB" ] || { echo "PREFLIGHT FAIL: no antiSMASH db at $DB" >&2; exit 1; }
[ -d "$MODELS/TXSScan" ] || { echo "PREFLIGHT FAIL: no TXSScan models" >&2; exit 1; }
echo "  6 Set B genomes, $(grep -c '^>' "$IN/setB_pC3_all.faa") proteins; all inputs ok"

source "$M/etc/profile.d/conda.sh"

# ------------------------------------------------------------- InterProScan
# -dp disables the EBI precalculated-match lookup: slower, but reproducible from
# on-disk data alone. All 17 member databases run (-appl omitted).
echo "=== [1/3] InterProScan ($(date +%H:%M)) ==="
"$IPS/interproscan.sh" -i "$IN/setB_pC3_all.faa" -b "$OUT/setB_pC3_ips" \
    -f TSV,GFF3 -goterms -pa -iprlookup -dp -cpu "$T" -T "$TMPDIR"
echo "  InterProScan rows: $(wc -l < "$OUT/setB_pC3_ips.tsv")"

# ---------------------------------------------------------------- antiSMASH
echo "=== [2/3] antiSMASH ($(date +%H:%M)) ==="
conda activate "$W/envs/antismash7"
for g in $MEM; do
    echo "  --- $g ---"
    antismash --taxon bacteria --genefinding-tool none --cpus "$T" \
        --databases "$DB" --cb-general --cb-knownclusters --asf --pfam2go \
        --output-dir "$OUT/antismash_$g" "$IN/${g}_pC3.gbk"
done
conda deactivate

# -------------------------------------------------------------- MacSyFinder
echo "=== [3/3] MacSyFinder ($(date +%H:%M)) ==="
conda activate "$W/envs/macsy"
for g in $MEM; do
    echo "  --- $g ---"
    macsyfinder --db-type ordered_replicon --sequence-db "$IN/${g}_pC3.faa" \
        --models TXSScan/bacteria/diderm all --models-dir "$MODELS" \
        -o "$OUT/macsy_$g" -w "$T"
done
conda deactivate

echo "=== ALL THREE SCREENS DONE ($(date +%H:%M)) ==="
ls "$OUT"
exit
