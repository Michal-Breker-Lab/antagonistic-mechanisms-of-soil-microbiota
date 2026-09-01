#!/bin/bash
# Search the 142 newly-called genomes for the 10 toxin/immunity loci, so every
# tip of the bac120 tree has a measured value instead of a "never assessed" cell.
#
# EVERYTHING that could make these results incomparable with the existing table
# is copied from s27_search_toxins.py rather than re-chosen:
#   * the same 12 query sequences (queries_toxins12.faa) = 10 LOCI, two of which
#     are searched at two scales (full-length + C-terminal tip). Counting 12 as
#     12 independent proteins double-counts MF6_003684 and MF6_004284.
#   * the same -outfmt field list, in the same order
#   * -evalue 1e-5, -max_target_seqs 1000000, -comp_based_stats 2
#   * the same tier rule: tier1 = pident >= 60 AND qcovs >= 70
#
# MULTI-DATABASE SYNTAX IS LOAD-BEARING, and it is the subtle one. BLAST sums the
# effective database sizes across a quoted -db list, so E-values computed here
# stay on the same scale as s27's. Searching the 142 new genomes as their own
# small database would inflate their E-values relative to the 771 and make the
# 1e-5 gate mean something different for them than for everyone else.
set -euo pipefail
C=/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3
W=/mnt/LargeStorageNoBackup/Moshea/burk_bac120
BIN=$C/envs/bakta/bin
export PATH="$BIN:$PATH"
export TMPDIR=/mnt/LargeStorageNoBackup/Moshea/tmp
mkdir -p "$TMPDIR" "$W/search_toxins" "$W/logs"
cd "$W"
CPUS=${CPUS:-48}

FMT="6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qcovs slen qlen"

echo "=== preflight ==="
blastp -version | head -1
[ -s "$C/queries_toxins12.faa" ] || { echo "FAIL: no query file" >&2; exit 1; }
echo "  queries: $(grep -c '^>' "$C/queries_toxins12.faa")"

# --- build a database of the NEW genomes only --------------------------------
NEW=$W/search_toxins/new_pyrodigal.faa
if [ ! -s "$NEW.pin" ] && [ ! -s "$NEW.pdb" ]; then
    : > "$NEW"
    n=0
    while read -r a; do
        cat "$C/pyrodigal_prot/$a.faa" >> "$NEW"; n=$((n+1))
    done < need_prot.txt
    echo "  new-genome FASTA: $n genomes, $(grep -c '^>' "$NEW") proteins"
    makeblastdb -in "$NEW" -dbtype prot -out "$NEW" -title burk_new >/dev/null
fi

# --- search: new DB + the two originals, so E-values stay comparable ---------
TSV=$W/search_toxins/blastp_hits_new.tsv
if [ ! -s "$TSV" ]; then
    blastp -query "$C/queries_toxins12.faa" \
           -db "$NEW $C/search_rhs/all_pyrodigal.faa $C/pyrodigal_mf6/MF6.faa" \
           -outfmt "$FMT" -evalue 1e-5 -max_target_seqs 1000000 \
           -num_threads "$CPUS" -comp_based_stats 2 -out "$TSV"
fi
echo "  blastp rows: $(wc -l < "$TSV")"

echo "=== tier assignment ==="
python3 - <<'PY'
import csv, os
from collections import defaultdict
W = "/mnt/LargeStorageNoBackup/Moshea/burk_bac120"
COLS = ("qseqid sseqid pident length nident mismatch gapopen qstart qend sstart "
        "send evalue bitscore qcovhsp qcovs slen qlen").split()
QUERIES = ["MF6_001079", "MF6_002734", "MF6_003184", "MF6_003684_2514_2867",
           "MF6_003684_full", "MF6_003686", "MF6_003843",
           "MF6_004284_1434_1538", "MF6_004284_full", "MF6_004285",
           "MF6_004947", "MF6_006318"]
# Figure 11/12 draw these ten; the two _Cterm names are the C-terminal queries.
RING = {"MF6_001079": "MF6_001079", "MF6_002734": "MF6_002734",
        "MF6_003184": "MF6_003184", "MF6_003684_2514_2867": "MF6_003684_Cterm",
        "MF6_003686": "MF6_003686", "MF6_003843": "MF6_003843",
        "MF6_004284_1434_1538": "MF6_004284_Cterm", "MF6_004285": "MF6_004285",
        "MF6_004947": "MF6_004947", "MF6_006318": "MF6_006318"}


def tier(pid, qcov):
    if pid >= 60 and qcov >= 70:
        return 1
    if pid >= 60 and qcov >= 30:
        return 2
    if 40 <= pid < 60 and qcov >= 70:
        return 3
    return 4


best = defaultdict(lambda: 9)          # (acc, ring_locus) -> best tier
bad = 0
with open(f"{W}/search_toxins/blastp_hits_new.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) != len(COLS):
            bad += 1
            continue
        rec = dict(zip(COLS, f))
        parts = rec["sseqid"].split("|")
        if len(parts) != 6:
            bad += 1
            continue
        q = rec["qseqid"]
        if q not in RING:
            continue
        acc = parts[0]
        t = tier(float(rec["pident"]), float(rec["qcovs"]))
        k = (acc, RING[q])
        if t < best[k]:
            best[k] = t
if bad:
    raise SystemExit(f"ERROR: {bad} malformed rows - do not interpret these results")

searched = [l.strip() for l in open(f"{W}/all_bac120.txt") if l.strip()]
out = f"{W}/toxin12_carriers_bac120.tsv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "locus", "best_tier", "tier1_carrier"])
    for acc in searched:
        for locus in sorted(set(RING.values())):
            t = best.get((acc, locus), 0)
            w.writerow([acc, locus, t if t != 9 else "", str(t == 1)])
n1 = sum(1 for k, v in best.items() if v == 1)
print(f"  genomes searched: {len(searched)}")
print(f"  (genome, locus) tier-1 calls: {n1}")
print(f"  wrote {out}")
PY
echo "=== TOXIN SEARCH OK ==="
