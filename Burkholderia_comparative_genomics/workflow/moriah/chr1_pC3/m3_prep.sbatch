#!/bin/bash
#SBATCH --job-name=c3_prep
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/m3_prep.%j.out
#SBATCH --error=/sci/backup/ofinkel/moshea/burkholderia_c3/logs/m3_prep.%j.err
# NO --partition: the -k variants are account-restricted; SLURM picks via ofinkel.
# 32 cores x 128 GB = 4 GB/core -> base billing tier ($0.008/core-h).

set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
THREADS=${SLURM_CPUS_PER_TASK:-32}

# Node /tmp is per-user quota'd and blew up a bakta array before; keep all
# scratch on the lab filesystem and clean it up on any exit path.
mkdir -p "$W/tmp"
export TMPDIR
TMPDIR=$(mktemp -d "$W/tmp/prep.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

source "$M/etc/profile.d/conda.sh"

SKANI="$M/envs/gtdbtk/bin/skani"
PY="$M/bin/python"     # the pyrodigal launcher has a stale /sci/home shebang;
                       # the module itself imports fine under base python.

echo "=== preflight (RUN the tools, do not just test -x) ==="
"$SKANI" --version
"$PY" -c "import pyrodigal, sys; print('pyrodigal', pyrodigal.__version__)"
"$M/bin/makeblastdb" -version | head -1
echo "TMPDIR=$TMPDIR  THREADS=$THREADS"

mkdir -p "$W/results" "$W/proteins"

# ---------------------------------------------------------------- Stage A: ANI
echo; echo "=== Stage A: skani all-vs-all ==="
LIST="$TMPDIR/fasta_paths.txt"
: > "$LIST"
while read -r acc; do
  [ -z "$acc" ] && continue
  case "$acc" in
    MF6) echo "$W/MF6.fna" >> "$LIST" ;;
    MF7) echo "$W/MF7.fna" >> "$LIST" ;;
    *)   echo "$W/genomes/$acc.fna" >> "$LIST" ;;
  esac
done < "$W/genome_list_full.txt"
echo "genomes for skani: $(wc -l < "$LIST")"

# -E gives a row-per-pair edge list (same shape as `skani dist`), which is what
# the 99% clone clustering consumes.  Default --min-af 15 is kept so distant
# cross-genus pairs drop out instead of bloating the file.
"$SKANI" triangle -l "$LIST" -t "$THREADS" -E --short-header \
    -o "$W/results/ani_sparse.tsv"
echo "skani rows: $(wc -l < "$W/results/ani_sparse.tsv")"

# ---------------------------------------------------- Stage B: pyrodigal + DB
echo; echo "=== Stage B: pyrodigal proteins ==="
"$PY" "$W/s23_pyrodigal.py" \
    --genomes "$W/genomes" --list "$W/genome_list_full.txt" \
    --mf6 "$W/MF6.fna" --mf7 "$W/MF7.fna" \
    --outdir "$W/proteins" --jobs "$THREADS"

echo "faa files: $(ls "$W/proteins"/*.faa | wc -l)"
cat "$W/proteins"/*.faa > "$W/results/all_proteins.faa"
echo "total proteins: $(grep -c '^>' "$W/results/all_proteins.faa")"

"$M/bin/makeblastdb" -in "$W/results/all_proteins.faa" -dbtype prot \
    -out "$W/results/blastdb/all_proteins" -title all_proteins

echo; echo "=== DONE ==="
exit 0
