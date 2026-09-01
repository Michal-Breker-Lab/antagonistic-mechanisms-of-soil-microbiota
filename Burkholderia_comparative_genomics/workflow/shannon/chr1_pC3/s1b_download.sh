#!/usr/bin/env bash
# Stage 1b - download genome FASTA for every distinct assembly.
# Batched at 100 accessions so a single network failure costs one batch, not
# the whole 6 GB. Re-runnable: batches whose .done marker exists are skipped.
set -uo pipefail
W=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
export TMPDIR=$W/tmp
DATASETS=$W/envs/burk/bin/datasets
mkdir -p $W/raw/zips $W/genomes

split -l 100 -d -a 3 $W/metadata/accessions.txt $W/raw/batch_
n=0
for b in $W/raw/batch_*; do
  id=$(basename "$b")
  [ -f "$W/raw/zips/$id.done" ] && { echo "skip $id"; continue; }
  $DATASETS download genome accession --inputfile "$b" \
      --include genome --no-progressbar \
      --filename "$W/raw/zips/$id.zip" > "$W/logs/dl_$id.log" 2>&1
  if [ $? -eq 0 ] && [ -s "$W/raw/zips/$id.zip" ]; then
    unzip -o -q "$W/raw/zips/$id.zip" -d "$W/raw/zips/$id" \
      && touch "$W/raw/zips/$id.done"
    n=$((n+1)); echo "ok $id"
  else
    echo "FAIL $id"
  fi
done
echo "batches downloaded this run: $n"

# Flatten to one <accession>.fna per genome with the accession as filename.
echo "=== flattening ==="
found=0
for d in $W/raw/zips/batch_*/; do
  [ -d "$d/ncbi_dataset/data" ] || continue
  for gd in "$d"ncbi_dataset/data/GC*/; do
    acc=$(basename "$gd")
    f=$(ls "$gd"*.fna 2>/dev/null | head -1)
    [ -n "$f" ] || continue
    [ -f "$W/genomes/$acc.fna" ] || cp "$f" "$W/genomes/$acc.fna"
    found=$((found+1))
  done
done
echo "genomes present: $(ls $W/genomes/*.fna 2>/dev/null | wc -l) (expected $(wc -l < $W/metadata/accessions.txt))"
du -sh $W/genomes
echo "DOWNLOAD_DONE"
