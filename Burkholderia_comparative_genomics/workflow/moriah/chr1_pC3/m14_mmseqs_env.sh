#!/bin/bash
# Build a minimal MMseqs2 env for the Set B clade-core clustering.
# Run via srun on a compute node -- never on the login/gateway node
# (typed memory: moriah-build-envs-on-compute-node).
set -euo pipefail
M=/sci/backup/ofinkel/moshea/miniforge3
ENVDIR=/sci/backup/ofinkel/moshea/burkholderia_c3/envs/mmseqs2
source "$M/etc/profile.d/conda.sh"
export TMPDIR=/sci/backup/ofinkel/moshea/burkholderia_c3/tmp
mkdir -p "$TMPDIR"

"$M/bin/mamba" create -y -p "$ENVDIR" -c conda-forge -c bioconda mmseqs2

# Verify on the filesystem + a real invocation, not on the installer's exit code
# (job 45911107 exited 120 on a cosmetic check after a successful install).
[ -x "$ENVDIR/bin/mmseqs" ] || { echo "FAIL: no mmseqs binary" >&2; exit 1; }
"$ENVDIR/bin/mmseqs" version
echo "MMSEQS_ENV_OK $ENVDIR"
exit
