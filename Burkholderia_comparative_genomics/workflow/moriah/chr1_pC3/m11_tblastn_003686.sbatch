#!/bin/bash
#SBATCH --job-name=c3_tbn
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=logs/m11_tbn.%j.out
#SBATCH --error=logs/m11_tbn.%j.err
# NO --partition (the -k variants are account-restricted).
#
# tblastn for the ONE locus a protein-level search cannot answer: MF6_003686 is
# a 120-aa ORF that Bakta does not call and pyrodigal truncates to 78 aa via an
# internal start, so MF6 scores 65% coverage against its own gene -- below the
# tier-1 floor. Searching six-frame translations of the assemblies instead.
# Flags reproduce TOOLS.md exactly except -num_threads (4 here, 48 originally;
# thread count cannot change a hit).
set -euo pipefail

W=/sci/backup/ofinkel/moshea/burkholderia_c3
M=/sci/backup/ofinkel/moshea/miniforge3
THREADS=${SLURM_CPUS_PER_TASK:-4}
cd "$W"

export TMPDIR=$(mktemp -d "$W/tmp/tbn.${SLURM_JOB_ID:-x}.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== preflight: CHECK THE INPUTS, then run the tools ==="
# BLAST+ lives in the base miniforge3, same as m4_search.sbatch used for blastp.
# There is no `blast` conda env on Moriah -- call the binaries by absolute path.
TBLASTN="$M/bin/tblastn"; MAKEDB="$M/bin/makeblastdb"
[ -x "$TBLASTN" ] || { echo "PREFLIGHT FAIL: $TBLASTN missing" >&2; exit 1; }
[ -x "$MAKEDB" ]   || { echo "PREFLIGHT FAIL: $MAKEDB missing" >&2; exit 1; }
"$TBLASTN" -version | head -1
[ -s "$W/queries_toxin12.faa" ] || { echo "PREFLIGHT FAIL: missing queries_toxin12.faa" >&2; exit 1; }
NG=$(ls "$W"/genomes/*.fna 2>/dev/null | wc -l)
[ "$NG" -eq 771 ] || { echo "PREFLIGHT FAIL: $NG genomes in genomes/, expected 771" >&2; exit 1; }
for f in MF6.fna MF7.fna; do
    [ -s "$W/$f" ] || { echo "PREFLIGHT FAIL: missing $f" >&2; exit 1; }
done
echo "  inputs ok: 771 NCBI genomes + MF6 + MF7"

# --- extract the single query ----------------------------------------------
Q="$TMPDIR/q_003686.faa"
"$M/bin/python" - "$W/queries_toxin12.faa" "$Q" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
keep, out = False, []
for line in open(src):
    if line.startswith(">"):
        keep = line[1:].split()[0] == "MF6_003686"
    if keep:
        out.append(line)
if not out:
    sys.exit("FAIL: MF6_003686 not found in the query FASTA")
open(dst, "w").writelines(out)
aa = sum(len(l.strip()) for l in out[1:])
print(f"  query MF6_003686: {aa} aa")
PY

# --- nucleotide DB, headers as <accession>|<contig> -------------------------
DB="$TMPDIR/all_genomes.fna"
"$M/bin/python" - "$W" "$DB" <<'PY'
import sys, glob, os
W, dst = sys.argv[1], sys.argv[2]
n_seq = n_gen = 0
with open(dst, "w") as out:
    files = sorted(glob.glob(os.path.join(W, "genomes", "*.fna")))
    files += [os.path.join(W, "MF6.fna"), os.path.join(W, "MF7.fna")]
    for p in files:
        acc = os.path.basename(p)[:-4]
        n_gen += 1
        for line in open(p):
            if line.startswith(">"):
                out.write(f">{acc}|{line[1:].split()[0]}\n")
                n_seq += 1
            else:
                out.write(line)
print(f"  db: {n_gen} genomes, {n_seq} contigs")
PY
"$MAKEDB" -in "$DB" -dbtype nucl -out "$TMPDIR/allg" >/dev/null

# --- the search -------------------------------------------------------------
COLS="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qlen"
"$TBLASTN" -query "$Q" -db "$TMPDIR/allg" -db_gencode 11 -evalue 1e-5 \
        -max_target_seqs 1000000 -num_threads "$THREADS" \
        -outfmt "$COLS" -out "$W/results/toxin12_tblastn_003686.tsv"

echo "=== done: $(wc -l < "$W/results/toxin12_tblastn_003686.tsv") hits ==="
head -3 "$W/results/toxin12_tblastn_003686.tsv"
exit
