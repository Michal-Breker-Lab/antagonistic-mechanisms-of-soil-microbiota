#!/bin/bash
#SBATCH --job-name=burk_dl
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/download_%j.out
#SBATCH --error=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/download_%j.err
# No --partition: the -k variants are account-restricted and SLURM picks by account.
set -uo pipefail
W=/sci/backup/ofinkel/moshea/burkholderia_c3
source /sci/backup/ofinkel/moshea/miniforge3/etc/profile.d/conda.sh
conda activate datasets

# Node /tmp is per-user quota'd; a previous bakta array lost 38% of its tasks to
# EDQUOT because of it. Everything temporary goes on the lab filesystem.
export TMPDIR=$(mktemp -d "$W/tmp/dl.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

mkdir -p "$W/genomes" "$W/raw"
grep -vE '^(MF6|MF7)$' "$W/genome_list_full.txt" > "$TMPDIR/accs.txt"
echo "accessions to fetch: $(wc -l < "$TMPDIR/accs.txt")"

# Batch at 100: a single HTTP/2 stream failure then costs one batch, not 6 GB.
split -l 100 -d -a 3 "$TMPDIR/accs.txt" "$TMPDIR/batch_"
for b in "$TMPDIR"/batch_*; do
  id=$(basename "$b")
  [ -f "$W/raw/$id.done" ] && { echo "skip $id"; continue; }
  datasets download genome accession --inputfile "$b" \
      --include genome --no-progressbar --filename "$W/raw/$id.zip" \
      > "$W/logs/dl_$id.log" 2>&1
  rc=$?
  # NCBI truncates streams while still exiting 0 -- the log is the real check.
  if grep -qi "INTERNAL_ERROR\|stream error" "$W/logs/dl_$id.log"; then
    echo "STREAM ERROR on $id - not marking done"; continue
  fi
  if [ $rc -eq 0 ] && [ -s "$W/raw/$id.zip" ]; then
    unzip -o -q "$W/raw/$id.zip" -d "$W/raw/$id" && touch "$W/raw/$id.done"
    echo "ok $id"
  else
    echo "FAIL $id (rc=$rc)"
  fi
done

echo "=== flattening to <accession>.fna ==="
for gd in "$W"/raw/batch_*/ncbi_dataset/data/GC*/; do
  [ -d "$gd" ] || continue
  acc=$(basename "$gd")
  f=$(ls "$gd"*.fna 2>/dev/null | head -1)
  [ -n "$f" ] && [ ! -s "$W/genomes/$acc.fna" ] && cp "$f" "$W/genomes/$acc.fna"
done

n=$(ls "$W"/genomes/*.fna 2>/dev/null | wc -l)
want=$(wc -l < "$TMPDIR/accs.txt")
echo "genomes present: $n / $want"
missing=0
while read -r a; do [ -s "$W/genomes/$a.fna" ] || { echo "MISSING $a"; missing=$((missing+1)); }; done < "$TMPDIR/accs.txt"
echo "missing: $missing"
du -sh "$W/genomes"
[ "$missing" -eq 0 ] && echo "DOWNLOAD_OK" || echo "DOWNLOAD_INCOMPLETE"
exit 0
