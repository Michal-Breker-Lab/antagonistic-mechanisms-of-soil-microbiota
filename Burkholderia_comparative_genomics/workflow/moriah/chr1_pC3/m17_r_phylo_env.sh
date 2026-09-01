#!/bin/bash
# R phylogenetics env for the section 7 comparative analyses (stochastic character
# mapping, Fritz & Purvis D, phylogenetic logistic regression). Moriah's base R
# is broken ("Rscript execution error"); the anvio-9 env's R 4.5.2 works but
# carries none of these packages, so build a dedicated env.
# Run via srun on a compute node, never the login node.
set -euo pipefail
M=/sci/backup/ofinkel/moshea/miniforge3
ENVDIR=/sci/backup/ofinkel/moshea/burkholderia_c3/envs/rphylo
export TMPDIR=/sci/backup/ofinkel/moshea/burkholderia_c3/tmp
mkdir -p "$TMPDIR"

"$M/bin/mamba" create -y -p "$ENVDIR" -c conda-forge \
    r-base r-ape r-phytools r-caper r-geiger r-phylolm r-nlme

# Verify by loading, not by trusting the installer's exit code.
"$ENVDIR/bin/Rscript" -e '
for (p in c("ape","phytools","caper","geiger","phylolm")) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(sprintf("%-10s %s %s\n", p, ifelse(ok,"OK","FAIL"),
              ifelse(ok, as.character(packageVersion(p)), "")))
  if (!ok) quit(status=1)
}'
echo "RPHYLO_ENV_OK $ENVDIR"
exit
