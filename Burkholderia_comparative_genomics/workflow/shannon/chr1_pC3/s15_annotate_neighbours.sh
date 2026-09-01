#!/usr/bin/env bash
# Annotate the 3 MF6 near-neighbours that dereplication removed before Stage 6,
# using EXACTLY the Stage 6 flags so their GFF3s are interchangeable with the rest.
#
#   GCF_016899425.1  MS389            ANI 99.71 to MF6  <- MF6's closest relative
#   GCF_053209605.1  B. sola BC00106  ANI 96.11
#   GCF_000203955.1  B. cenocepacia HI2424  ANI 94.63
#
# Their c3 calls in c3_calls_all_genomes.tsv read "inherited" (propagated from
# their dereplication representative), never measured. Annotating them lets them
# enter the pangenome directly.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
DB=/mnt/LargeStorageNoBackup/Datasets/Moshea/Databases/bakta_db/db
# Bakta resolves tRNAscan-SE/aragorn/cmscan/pilercr/diamond from PATH, not from
# its own prefix -- absolute-path invocation alone fails with "tRNAscan-SE not found".
export PATH=$W/envs/bakta/bin:$PATH
BAKTA=$W/envs/bakta/bin/bakta
mkdir -p "$W/logs/bakta"

for A in GCF_016899425.1 GCF_053209605.1 GCF_000203955.1; do
  if [ -s "$W/annot/$A/$A.gff3" ]; then
    echo "SKIP $A (already annotated)"
    continue
  fi
  mkdir -p "$W/annot/$A"
  echo "=== bakta $A  $(date +%H:%M:%S) ==="
  # --keep-contig-headers is REQUIRED: the GFF3 subsetting step matches contigs
  # by their original NCBI accession and Bakta renames them without it.
  # --skip-plot saves several minutes producing circular figures we never use.
  if $BAKTA --db "$DB" --output "$W/annot/$A" --prefix "$A" \
            --genus Burkholderia --species sp. \
            --keep-contig-headers --skip-plot \
            --threads 24 --force "$W/genomes/$A.fna" \
            > "$W/logs/bakta/$A.log" 2>&1; then
    echo "OK $A  $(date +%H:%M:%S)"
  else
    echo "FAIL $A -- see $W/logs/bakta/$A.log"
    tail -5 "$W/logs/bakta/$A.log"
  fi
done
echo "S15_ANNOT_DONE"
exit
