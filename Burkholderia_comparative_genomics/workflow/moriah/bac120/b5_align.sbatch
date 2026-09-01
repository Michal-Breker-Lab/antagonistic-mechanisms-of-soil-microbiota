#!/bin/bash
#SBATCH --job-name=bac120_aln
#SBATCH --partition=glacier
#SBATCH --exclude=glacier-11
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=logs/b5_align.%j.out
#SBATCH --error=logs/b5_align.%j.err
#
# 128G over 32 cores = 4 GB/core -> base billing tier.
# --partition / --exclude: see the note in b3_download.sbatch.
set -euo pipefail
W=/sci/backup/ofinkel/moshea/burk_bac120
M=/sci/backup/ofinkel/moshea/miniforge3
cd "$W"
export TMPDIR=$(mktemp -d "$W/tmp/aln.${SLURM_JOB_ID:-x}.XXXXXX"); trap 'rm -rf "$TMPDIR"' EXIT

# Same runtime env resolution as b4 -- see the long comment there.
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
[ -s identify/identify/gtdbtk.bac120.markers_summary.tsv ] \
  || { echo "PREFLIGHT FAIL: identify output missing -- run b4 first" >&2; exit 1; }

echo "=== gtdbtk align (--skip_gtdb_refs) ==="
# --skip_gtdb_refs is LOAD-BEARING. Without it thousands of GTDB reference genomes
# join the MSA and the tree stops being about our 792 genomes.
gtdbtk align --identify_dir identify --out_dir align \
       --skip_gtdb_refs --cpus "${SLURM_CPUS_PER_TASK:-32}"

echo "=== verify ==="
A=align/align/gtdbtk.bac120.user_msa.fasta.gz
[ -s "$A" ] || { echo "FAIL: no user MSA written at $A" >&2; exit 1; }
gunzip -c "$A" > msa_bac120.faa
TAXA=$(grep -c '^>' msa_bac120.faa)
COLS=$(awk '/^>/{if(n){print length(s); exit} n=1; next} {s=s $0}' msa_bac120.faa)
echo "  taxa: $TAXA"
echo "  columns: $COLS"
[ "$TAXA" -gt 700 ] || { echo "FAIL: only $TAXA taxa in the MSA" >&2; exit 1; }
grep -q "^>MF7" msa_bac120.faa || { echo "FAIL: MF7 absent from the MSA" >&2; exit 1; }
grep -q "^>MF6" msa_bac120.faa || { echo "FAIL: MF6 absent from the MSA" >&2; exit 1; }
echo "  MF6 and MF7 both in the alignment"
# Taxa < 792 is NOT an error -- gtdbtk align drops genomes with too few markers.
# Report the shortfall explicitly; a silent drop reads as "everything passed".
if [ "$TAXA" -ne 792 ]; then
    echo "  NOTE: $(( 792 - TAXA )) genome(s) dropped by align for insufficient markers;"
    echo "        b6_qc_markers.py enumerates them. This is reported, not hidden."
fi
echo "=== ALIGN OK ==="
exit
