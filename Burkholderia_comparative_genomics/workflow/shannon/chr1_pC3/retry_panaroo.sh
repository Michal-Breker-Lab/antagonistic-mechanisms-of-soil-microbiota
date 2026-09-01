export TMPDIR=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/tmp; export CONDA_PKGS_DIRS=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/tmp/pkgs
~/miniconda/bin/mamba create -y -p /mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/envs/panaroo -c conda-forge -c bioconda panaroo 2>&1 | tail -8
echo '--- version ---'; /mnt/LargeStorageNoBackup/Moshea/burkholderia_c3/envs/panaroo/bin/panaroo --version 2>&1 | head -2
echo PANAROO_ENV_DONE
