#!/usr/bin/env bash
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
export CONDA_PKGS_DIRS=$W/tmp/pkgs      # keep tarballs off the 93%-full /
mkdir -p "$CONDA_PKGS_DIRS"
MAMBA=~/miniconda/bin/mamba

echo "=== env 1/3: burk (core toolkit) ==="
$MAMBA create -y -p $W/envs/burk -c conda-forge -c bioconda \
  python=3.11 ncbi-datasets-cli skani seqkit mmseqs2 hmmer mafft "iqtree>=2.2" \
  biopython pandas numpy scipy matplotlib 2>&1 | tail -5

echo "=== env 2/3: bakta ==="
$MAMBA create -y -p $W/envs/bakta -c conda-forge -c bioconda bakta 2>&1 | tail -5

echo "=== env 3/3: panaroo ==="
$MAMBA create -y -p $W/envs/panaroo -c conda-forge -c bioconda panaroo 2>&1 | tail -5

echo "=== VERSIONS ==="
$W/envs/burk/bin/datasets --version 2>&1 | head -1
$W/envs/burk/bin/skani --version 2>&1 | head -1
$W/envs/burk/bin/seqkit version 2>&1 | head -1
$W/envs/burk/bin/mmseqs version 2>&1 | head -1
$W/envs/burk/bin/hmmsearch -h 2>&1 | grep -m1 HMMER
$W/envs/burk/bin/mafft --version 2>&1 | head -1
$W/envs/burk/bin/iqtree2 --version 2>&1 | head -1 || $W/envs/burk/bin/iqtree --version 2>&1 | head -1
$W/envs/bakta/bin/bakta --version 2>&1 | head -1
$W/envs/panaroo/bin/panaroo --version 2>&1 | head -1
df -h / /mnt/LargeStorageNoBackup | tail -2
echo "SETUP_DONE"
