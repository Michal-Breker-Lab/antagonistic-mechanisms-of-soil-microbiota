#!/usr/bin/env bash
# Core gene families of MF6's closest relatives, clustered on EXACTLY those genomes.
#
# Set A = the 10 closest by chromosome-1 ANI + MF6 (11 genomes). Spans the ANI
#         break at rank 4, so it mixes B. sola with assorted B. cenocepacia.
# Set B = the B. sola clade (ANI >= 95%) + MF6 (5 genomes).
#
# Clustering on the subset (rather than subsetting the 140/306-genome Rtab) is the
# direct answer to "what is the core of these genomes". The subset-of-the-big-Rtab
# version is computed separately by s16b_neighbour_core_from_rtab.py and reported
# alongside, because only that one is comparable to the main pangenome's families.
#
# NOTE: PPanGGOLiN's statistical persistent/shell/cloud partitioning is not
# reliable at n=5-11 genomes; for these sets report the PREVALENCE-THRESHOLD core
# from the Rtab, not the persistent partition.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
PPGENV=/mnt/scratch/within_genus_hgt/env_ppanggolin
PPG=$PPGENV/bin/ppanggolin
export PATH=$PPGENV/bin:$PATH          # ppanggolin resolves mmseqs from PATH
CPU=${CPU:-16}

SET_A="GCF_016899425.1 GCF_905400185.1 GCF_053038975.1 GCF_053209605.1 \
GCF_003966315.1 GCF_001718515.1 GCF_014211915.1 GCF_003076415.1 \
GCF_000019505.1 GCF_000203955.1 MF6"
SET_B="GCF_016899425.1 GCF_905400185.1 GCF_053038975.1 GCF_053209605.1 MF6"

mkdir -p "$W/ppanggolin/neighbours" "$W/logs"

# build an --anno list, drawing each genome from the main GFF3 dir or the _plus
# dir (the 3 genomes annotated late). Genomes with no c3 simply have no gff_c3
# file and are silently absent from the c3 sets -- that is correct, not an error.
build_list () {                 # $1 = setname  $2 = replicon (c3|chr1)  $3.. = accessions
  local SETNAME=$1 REP=$2; shift 2
  local L=$W/ppanggolin/neighbours/list_${SETNAME}_${REP}.tsv
  : > "$L"
  for a in "$@"; do
    for d in "$W/pangenome/gff_$REP" "$W/ppanggolin/gff_${REP}_plus"; do
      if [ -s "$d/$a.gff3" ]; then printf '%s\t%s\n' "$a" "$d/$a.gff3" >> "$L"; break; fi
    done
  done
  echo "$L: $(wc -l < "$L") genomes"
}

run_set () {                    # $1 = setname  $2 = replicon  $3 = identity
  local SETNAME=$1 REP=$2 ID=$3
  local L=$W/ppanggolin/neighbours/list_${SETNAME}_${REP}.tsv
  local O=$W/ppanggolin/neighbours/${SETNAME}_${REP}_id${ID/./}
  local H=$O/pangenome.h5
  mkdir -p "$O"
  echo "=== START $SETNAME $REP id=$ID ($(wc -l < "$L") genomes) ==="
  $PPG annotate --anno "$L" -o "$O" -f -c "$CPU"                   || return 1
  $PPG cluster -p "$H" --identity "$ID" --coverage 0.8 --mode 1 -c "$CPU" || return 1
  $PPG graph -p "$H"                                               || return 1
  # partitioning is unreliable at this n; allow it to fail without killing the run
  $PPG partition -p "$H" -c "$CPU" || echo "WARN partition failed (expected at small n)"
  $PPG write_pangenome -p "$H" -o "$O/out" -f --Rtab --csv --stats || return 1
  echo "=== DONE $SETNAME $REP id=$ID ==="
}

for REP in c3 chr1; do
  build_list setA "$REP" $SET_A
  build_list setB "$REP" $SET_B
done

for SETNAME in setA setB; do
  for REP in c3 chr1; do
    for ID in 0.60 0.80; do
      run_set "$SETNAME" "$REP" "$ID" 2>&1 \
        | tee -a "$W/logs/s16_${SETNAME}_${REP}_${ID/./}.log"
    done
  done
done
echo "S16_ALL_DONE"
exit
