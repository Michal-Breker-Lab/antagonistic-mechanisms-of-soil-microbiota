#!/bin/bash
# tblastn rescue for MF6_003686 on the newly-added genomes.
#
# WHY THIS LOCUS ONLY, and why it cannot be skipped: MF6_003686 is a 120 aa
# query, and blastp against the pyrodigal proteins caps its coverage at qcovs=65
# -- below the tier-1 bar of 70 -- because pyrodigal's ORF call for it is
# truncated. Run blastp alone and MF6 fails to carry its OWN locus, which is how
# this was caught: a naive re-run disagreed with the retained table on exactly 9
# (genome, locus) cells, all of them MF6_003686, MF6 and MF7 among them.
# The retained pipeline handled it with a dedicated tblastn against NUCLEOTIDE,
# which recovers the full 120 aa at 100% identity (tables/toxin12_tblastn_003686.tsv).
# The carrier rule is unchanged -- tier1 = pident >= 60 AND qcov >= 70 -- only the
# search program differs. Verified against the retained table: carriers there span
# pident 89.9-100 / qcov 76-100, non-carriers fail at pident 58.9 or qcov 63.
#
# The other nine loci need no rescue: blastp reproduces the retained table on
# 6,541 of 6,550 comparable cells (99.86%), and all nine misses were this locus.
set -euo pipefail
C=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
BIN=$C/envs/bakta/bin
export PATH="$BIN:$PATH"
export TMPDIR=/mnt/LargeStorageNoBackup/Moshea/tmp
mkdir -p "$TMPDIR" "$W/search_toxins" "$W/logs"
cd "$W"
CPUS=${CPUS:-48}

echo "=== preflight ==="
tblastn -version | head -1
Q=$W/search_toxins/q_003686.faa
python3 - <<'PY'
import os
C = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
W = "/mnt/LargeStorageNoBackup/Moshea/burk_bac120"
keep, out, w = "MF6_003686", [], False
for line in open(f"{C}/queries_toxins12.faa"):
    if line.startswith(">"):
        w = line[1:].split()[0] == keep
    if w:
        out.append(line)
assert out, "MF6_003686 not found in queries_toxins12.faa"
os.makedirs(f"{W}/search_toxins", exist_ok=True)
open(f"{W}/search_toxins/q_003686.faa", "w").writelines(out)
print(f"  query: {sum(len(l.strip()) for l in out if not l.startswith('>'))} aa")
PY

# nucleotide DB of the NEW genomes, headers "<acc>|<contig>" to match the
# retained table's sseqid convention
NDB=$W/search_toxins/new_nucl.fna
if [ ! -s "$NDB.nin" ] && [ ! -s "$NDB.ndb" ]; then
    : > "$NDB"
    while read -r a; do
        awk -v acc="$a" '/^>/{print ">" acc "|" substr($1,2); next} {print}' \
            "$W/genomes/$a.fna" >> "$NDB"
    done < need_prot.txt
    makeblastdb -in "$NDB" -dbtype nucl -out "$NDB" -title burk_new_nucl >/dev/null
    echo "  nucleotide DB: $(grep -c '^>' "$NDB") contigs"
fi

FMT="6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovs qlen"
TSV=$W/search_toxins/tblastn_003686_new.tsv
if [ ! -s "$TSV" ]; then
    tblastn -query "$Q" -db "$NDB" -outfmt "$FMT" -evalue 1e-5 \
            -max_target_seqs 1000000 -num_threads "$CPUS" -out "$TSV"
fi
echo "  tblastn rows: $(wc -l < "$TSV")"

python3 - <<'PY'
import csv
W = "/mnt/LargeStorageNoBackup/Moshea/burk_bac120"
H = ("qseqid sseqid pident length qstart qend sstart send evalue bitscore "
     "qcovs qlen").split()
best = {}
for line in open(f"{W}/search_toxins/tblastn_003686_new.tsv"):
    f = line.rstrip("\n").split("\t")
    if len(f) != len(H):
        continue
    r = dict(zip(H, f))
    acc = r["sseqid"].split("|")[0]
    v = (float(r["pident"]), float(r["qcovs"]))
    if acc not in best or v > best[acc]:
        best[acc] = v
targets = [l.strip() for l in open(f"{W}/need_prot.txt") if l.strip()]
out = f"{W}/tblastn_003686_new_calls.tsv"
n = 0
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "locus", "best_pident", "best_qcov", "tier1_carrier"])
    for a in targets:
        p, q = best.get(a, (0.0, 0.0))
        t1 = (p >= 60 and q >= 70)
        n += t1
        w.writerow([a, "MF6_003686", f"{p:.3f}", f"{q:.1f}", str(t1)])
print(f"  genomes with a hit: {len(best)}   NEW tier-1 carriers: {n}")
print(f"  wrote {out}")
PY
echo "=== TBLASTN OK ==="
