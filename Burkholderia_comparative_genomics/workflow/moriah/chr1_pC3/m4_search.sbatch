#!/bin/bash
#SBATCH --job-name=c3_search
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/m4_search.%j.out
#SBATCH --error=logs/m4_search.%j.err
# NO --partition: the -k variants are account-restricted; SLURM picks by account.
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
PY="$M/bin/python"
THREADS=${SLURM_CPUS_PER_TASK:-32}
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/search.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight (RUN the tools, do not just test -x) ==="
"$M/bin/blastp" -version | head -1
"$PY" -c "import pyrodigal; print('pyrodigal', pyrodigal.__version__)"
echo "TMPDIR=$TMPDIR  THREADS=$THREADS"

# Inputs, not just tools. The previous run passed preflight and then died in the
# join step because results/contig_ranks.tsv was never staged to Moriah. Check
# every file the job reads BEFORE spending any core-hours on it.
for f in "$W/genome_list_full.txt" "$W/queries_toxin12.faa" "$W/queries_rhs.faa" \
         "$W/results/contig_ranks.tsv" "$W/MF6.fna" "$W/MF7.fna"; do
    [ -s "$f" ] || { echo "PREFLIGHT FAIL: missing or empty input $f" >&2; exit 1; }
    echo "  input ok: $f ($(wc -l < "$f") lines)"
done

# --- Stage A: re-call proteins WITH coordinates -----------------------------
# s23 previously wrote only >acc|contig|idx headers; the hit tables need
# start/end/strand. Re-run so proteins and coords come from one pass and cannot
# disagree. Deterministic (single mode, self-trained), so the .faa is unchanged.
if [ -s "$W/results/all_coords.tsv" ] && [ -s "$W/results/blastdb/all_proteins.pin" ] \
   && [ "${FORCE_REBUILD:-0}" != "1" ]; then
    echo "=== Stages A-C skipped: proteins, coords and blastdb already present ==="
    echo "  proteins: $(grep -c '>' "$W/results/all_proteins.faa")"
    echo "  coords:   $(( $(wc -l < "$W/results/all_coords.tsv") - 1 ))"
else
echo "=== Stage A: pyrodigal + coords ==="
"$PY" "$W/s23_pyrodigal.py" --genomes "$W/genomes" --list "$W/genome_list_full.txt" \
    --mf6 "$W/MF6.fna" --mf7 "$W/MF7.fna" --outdir "$W/proteins" --jobs "$THREADS"

echo "=== Stage B: concatenate ==="
cat "$W"/proteins/*.faa > "$W/results/all_proteins.faa"
head -1 "$(ls "$W"/proteins/*.coords.tsv | head -1)" > "$W/results/all_coords.tsv"
for f in "$W"/proteins/*.coords.tsv; do tail -n +2 "$f"; done >> "$W/results/all_coords.tsv"
echo "proteins: $(grep -c '>' "$W/results/all_proteins.faa")"
echo "coords:   $(( $(wc -l < "$W/results/all_coords.tsv") - 1 ))"

echo "=== Stage C: blast DB ==="
mkdir -p "$W/results/blastdb"
"$M/bin/makeblastdb" -in "$W/results/all_proteins.faa" -dbtype prot \
    -out "$W/results/blastdb/all_proteins" -title all_proteins

fi

# --- Stage D: searches ------------------------------------------------------
echo "=== Stage D: toxin12 (12 queries) ==="
"$PY" "$W/s24_search.py" --queries "$W/queries_toxin12.faa" \
    --db "$W/results/blastdb/all_proteins" --coords "$W/results/all_coords.tsv" \
    --ranks "$W/results/contig_ranks.tsv" --prefix toxin12 \
    --outdir "$W/results" --threads "$THREADS" --blastp "$M/bin/blastp" \
    ${REUSE_TOXIN12_RAW:+--raw "$W/results/toxin12_blastp_raw.tsv"}

echo "=== Stage E: RHS (3 queries) ==="
"$PY" "$W/s24_search.py" --queries "$W/queries_rhs.faa" \
    --db "$W/results/blastdb/all_proteins" --coords "$W/results/all_coords.tsv" \
    --ranks "$W/results/contig_ranks.tsv" --prefix rhs \
    --outdir "$W/results" --threads "$THREADS" --blastp "$M/bin/blastp"

echo "=== DONE ==="
exit 0
