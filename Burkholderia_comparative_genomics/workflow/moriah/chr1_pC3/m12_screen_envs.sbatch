#!/bin/bash
#SBATCH --job-name=c3_envs
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=0-08:00:00
#SBATCH --output=logs/m12_envs.%j.out
#SBATCH --error=logs/m12_envs.%j.err
# NO --partition (the -k variants are account-restricted).
#
# Build the two tool envs figs 7-8 need, on a COMPUTE node -- never the login /
# gateway node (typed memory: moriah-build-envs-on-compute-node).
#
# Shannon holds the originals and is down, so both are rebuilt here to the same
# versions TOOLS.md records:
#   antismash7  antiSMASH 7.1.0  + biopython <1.80 + the external binaries the
#               bioconda package does NOT pull in.
#   macsy       MacSyFinder 2.1.6 (2.1.4 CANNOT parse TXSScan 1.1.4 -- the models
#               declare grammar vers="2.1" and 2.1.4 only accepts "2.0")
# InterProScan needs no env: 5.76-107.0 with its 52 GB data/ is already at
# /sci/backup/ofinkel/moshea/Tools/interproscan-5.76-107.0/.
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
ENVS="$W/envs"
DB="$W/databases/antismash_db"
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/envs.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT
mkdir -p "$ENVS" "$DB"

source "$M/etc/profile.d/conda.sh"
echo "=== node $(hostname); mamba $($M/bin/mamba --version 2>&1 | head -1) ==="

# ---------------------------------------------------------------- antiSMASH
if [ ! -x "$ENVS/antismash7/bin/antismash" ]; then
    echo "=== building antismash7 ==="
    # biopython pinned <1.80: antiSMASH 7.1 calls Bio.Seq.UnknownSeq, removed in 1.80.
    "$M/bin/mamba" create -y -p "$ENVS/antismash7" -c bioconda -c conda-forge \
        "antismash=7.1.0" "biopython<1.80" python=3.10 \
        hmmer hmmer2 fasttree prodigal glimmerhmm diamond muscle blast
else
    echo "=== antismash7 already present ==="
fi
conda activate "$ENVS/antismash7"
antismash --version
python -c "import Bio; print('  biopython', Bio.__version__)"
for b in hmmscan hmmsearch fasttree prodigal glimmerhmm diamond muscle blastp; do
    command -v "$b" >/dev/null || { echo "PREFLIGHT FAIL: $b missing from antismash7" >&2; exit 1; }
done
echo "  all external binaries on PATH"

# antiSMASH databases (~9 GB). Idempotent: it skips what is already complete.
echo "=== antiSMASH databases -> $DB ==="
download-antismash-databases --database-dir "$DB"
du -sh "$DB"
conda deactivate

# ---------------------------------------------------------------- MacSyFinder
if [ ! -x "$ENVS/macsy/bin/macsyfinder" ]; then
    echo "=== building macsy ==="
    "$M/bin/mamba" create -y -p "$ENVS/macsy" -c bioconda -c conda-forge \
        "macsyfinder=2.1.6" hmmer
else
    echo "=== macsy already present ==="
fi
conda activate "$ENVS/macsy"
macsyfinder --version
V=$(macsyfinder --version 2>&1 | grep -oE "2\.[0-9]+\.[0-9]+" | head -1)
[ "$V" = "2.1.6" ] || { echo "PREFLIGHT FAIL: macsyfinder $V, need 2.1.6" >&2; exit 1; }

MODELS="$W/macsy_models"
mkdir -p "$MODELS"
if [ ! -d "$MODELS/TXSScan" ]; then
    echo "=== installing TXSScan models ==="
    macsydata install --target "$MODELS" TXSScan
fi
# Verify on the FILESYSTEM, not via a macsydata subcommand: `macsydata available`
# does not accept --models-dir and exits 120, which failed job 45911107 AFTER the
# install had already succeeded. (The binary is also deprecated for `msf_data`.)
DIDERM="$MODELS/TXSScan/definitions/bacteria/diderm"
[ -d "$DIDERM" ] || { echo "PREFLIGHT FAIL: $DIDERM missing" >&2; exit 1; }
for m in T1SS T2SS T3SS T4aP T6SSi pT4SSt; do
    [ -f "$DIDERM/$m.xml" ] || echo "  note: $m.xml not in this model set"
done
echo "  TXSScan diderm models: $(ls "$DIDERM"/*.xml | wc -l)"
conda deactivate

echo "=== ALL ENVS BUILT ==="
exit
