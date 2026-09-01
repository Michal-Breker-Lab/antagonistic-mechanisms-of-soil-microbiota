#!/bin/bash
# Build a GTDB-Tk env on Shannon that is VERSION-MATCHED to the r232 reference.
#
# WHY 2.7.x: gtdbtk hardcodes the mask filename per version. The reference on
# Shannon is release232 (masks/gtdb_r232_bac120.mask, metadata VERSION_DATA=r232),
# so a 2.4/2.6 toolkit would look for gtdb_r226_bac120.mask and die. That failure
# is LOUD (missing file), not silent -- which is why this pairing is safer than
# Moriah's r226 + gtdbtk 2.4.0, where the toolkit is one patch BELOW the minimum
# supported for its own reference.
#
# WHY --prefix ON LargeStorage: Shannon's / is 95% full (~23 GB). Conda envs
# default to $HOME on /, and a gtdbtk env plus deps is 1-2 GB. Do not risk it.
# Also install iqtree 3.1.1 here -- the build that made this project's two
# existing trees. Shannon's on-hand iqtree3 is 3.0.1; matching costs 2 minutes.
set -euo pipefail
M=/mnt/LargeStorageNoBackup/Datasets/Moshea/miniforge3
P=/mnt/LargeStorageNoBackup/Moshea/envs/gtdbtk-2.7.2
mkdir -p "$(dirname "$P")"
export MAMBA_ROOT_PREFIX="$M"
if [ -x "$P/bin/gtdbtk" ]; then echo "env already present at $P"; else
  "$M/bin/mamba" create -y --prefix "$P" -c conda-forge -c bioconda \
      "gtdbtk=2.7.2" "iqtree=3.1.1"
fi
echo "=== verify ==="
export GTDBTK_DATA_PATH=/mnt/scratch/gtdbtk_database/release232
"$P/bin/gtdbtk" --version 2>&1 | head -2
"$P/bin/iqtree" --version 2>&1 | head -1
for t in prodigal hmmsearch hmmalign; do
    printf "  %-10s %s\n" "$t" "$("$P/bin/$t" -h 2>&1 | head -1 | cut -c1-60)"
done
# Confirm the toolkit resolves the mask that actually exists.
# NOTE: gtdbtk derives the mask FILENAME from VERSION_DATA in the reference
# package's own metadata.txt -- it is not hardcoded per toolkit version. So an
# older toolkit would also find gtdb_r232_bac120.mask; pinning 2.7.2 buys a
# supported pairing, NOT protection from a missing-file crash. Assert anyway:
# a mask that fails to resolve would silently change which columns are kept.
"$P/bin/python" - <<'PY'
import os, sys, gtdbtk.config.common as m
cls = next(getattr(m, n) for n in dir(m)
           if isinstance(getattr(m, n), type) and hasattr(getattr(m, n), "MASK_BAC120"))
c = cls()
mask = os.path.join(c.MASK_DIR, c.MASK_BAC120)
print(f"  VERSION_DATA : {c.VERSION_DATA}")
print(f"  MASK_BAC120  : {mask}")
print(f"  mask exists  : {os.path.isfile(mask)}")
print(f"  pfam / tigr  : {os.path.isdir(c.PFAM_HMM_DIR)} {os.path.isfile(c.TIGRFAM_HMMS)}")
sys.exit(0 if os.path.isfile(mask) and os.path.isdir(c.PFAM_HMM_DIR) else 1)
PY
echo "=== ENV OK ==="
