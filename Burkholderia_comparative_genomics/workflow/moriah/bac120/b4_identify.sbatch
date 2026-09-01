#!/bin/bash
#SBATCH --job-name=bac120_id
#SBATCH --partition=glacier
#SBATCH --exclude=glacier-11
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=logs/b4_identify.%j.out
#SBATCH --error=logs/b4_identify.%j.err
#
# 256G over 64 cores = 4 GB/core, inside the <=8 GB/core base billing tier
# ($0.008/core-h). It is sized for the tier, not because identify needs the RAM.
# --partition / --exclude: see the note in b3_download.sbatch.
set -euo pipefail
W=/sci/backup/ofinkel/moshea/burk_bac120
M=/sci/backup/ofinkel/moshea/miniforge3
cd "$W"
export TMPDIR=$(mktemp -d "$W/tmp/id.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT

# --- resolve the GTDB-Tk env AT RUNTIME -------------------------------------
# GTDB r226 is officially paired with GTDB-Tk 2.4.1-2.6.1 and the pre-existing env
# is 2.4.0 -- one patch below the minimum. Rather than hard-code a guess at submit
# time, prefer a purpose-built 2.6.1 env if it exists and fall back to the older
# one, printing which was used. identify/align depend only on the marker HMMs,
# which are release-stable; it is classify/root/decorate (not used here) that read
# the taxonomy and RED files a mismatch would corrupt.
if [ -x "$M/envs/gtdbtk-2.6.1/bin/gtdbtk" ]; then
    ENV="$M/envs/gtdbtk-2.6.1"
elif [ -x "$M/envs/gtdbtk/bin/gtdbtk" ]; then
    ENV="$M/envs/gtdbtk"
else
    echo "PREFLIGHT FAIL: no gtdbtk env found under $M/envs" >&2; exit 1
fi
source "$M/etc/profile.d/conda.sh"
conda activate "$ENV"
export GTDBTK_DATA_PATH=/sci/backup/aerez/aerez/moshea/release226

echo "=== preflight ==="
echo "  env: $ENV"
gtdbtk --version | head -1
[ -d "$GTDBTK_DATA_PATH/markers" ] || { echo "PREFLIGHT FAIL: no markers/ in GTDBTK_DATA_PATH" >&2; exit 1; }
N=$(ls genomes/*.fna 2>/dev/null | wc -l)
[ "$N" -eq 792 ] || { echo "PREFLIGHT FAIL: $N genomes, expected 792" >&2; exit 1; }
# MF7 is the whole methodological point of this tree -- fail here, not 3 h in.
for g in MF6 MF7; do
    [ -s "genomes/$g.fna" ] || { echo "PREFLIGHT FAIL: genomes/$g.fna absent" >&2; exit 1; }
done
echo "  genomes: $N   TMPDIR: $TMPDIR   threads: ${SLURM_CPUS_PER_TASK:-64}"

echo "=== gtdbtk identify ==="
gtdbtk identify --genome_dir genomes --out_dir identify \
       --extension fna --cpus "${SLURM_CPUS_PER_TASK:-64}"

echo "=== verify ==="
S=identify/identify/gtdbtk.bac120.markers_summary.tsv
[ -s "$S" ] || { echo "FAIL: no bac120 marker summary written" >&2; exit 1; }
ROWS=$(( $(wc -l < "$S") - 1 ))
echo "  marker summary rows: $ROWS (expected 792)"
[ "$ROWS" -eq 792 ] || { echo "FAIL: $ROWS rows, expected 792" >&2; exit 1; }
grep -qP "^MF7\t" "$S" || { echo "FAIL: MF7 absent from the marker summary" >&2; exit 1; }
grep -qP "^MF6\t" "$S" || { echo "FAIL: MF6 absent from the marker summary" >&2; exit 1; }
awk -F'\t' 'NR==1 || $1=="MF6" || $1=="MF7"' "$S"
echo "=== IDENTIFY OK ==="
exit
