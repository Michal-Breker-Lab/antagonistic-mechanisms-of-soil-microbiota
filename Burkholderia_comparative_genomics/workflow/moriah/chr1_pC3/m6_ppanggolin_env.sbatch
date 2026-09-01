#!/bin/bash
#SBATCH --job-name=ppg_env
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/m6_ppg_env.%j.out
#SBATCH --error=logs/m6_ppg_env.%j.err
# NO --partition. Built on a COMPUTE node, never the login/gateway node.
set -euo pipefail

M=/sci/backup/ofinkel/moshea/miniforge3
ENV="$M/envs/ppanggolin"
export TMPDIR=$(mktemp -d "/sci/backup/ofinkel/moshea/burkholderia_c3/tmp/ppgenv.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

source "$M/etc/profile.d/conda.sh"
echo "=== building ppanggolin 2.2.6 (pinned to match the original run) ==="
"$M/bin/mamba" create -y -p "$ENV" -c conda-forge -c bioconda \
    "ppanggolin=2.2.6" "python=3.10"

echo "=== verify: RUN it, do not just test -x ==="
"$ENV/bin/ppanggolin" --version
"$ENV/bin/ppanggolin" cluster --help | head -20
echo "mmseqs bundled: $("$ENV/bin/mmseqs" version 2>/dev/null || echo 'not in env')"

echo "=== DONE ==="
exit 0
