#!/usr/bin/env python3
"""
s27 - search all 771 NCBI Burkholderia genomes + MF6 for 12 candidate toxin /
immunity queries supplied by Moshe (to_Omri_final_toxins.faa).

These 12 sequences are 10 LOCI: two of them are searched at two scales,
full-length and C-terminal tip, because that contrast is the whole point for
polymorphic toxins - a conserved delivery scaffold with a swappable warhead.
  MF6_004284_full / MF6_004284_1434_1538   CFFIHE_04281, Peptidase C39, pC3
  MF6_003684_full / MF6_003684_2514_2867   CFFIHE_03684, RHS repeat,     pC3
Counting them as 12 independent proteins in any statistic would double-count
those two loci.

Databases: the 771-genome database built in s24 PLUS a separate one-genome
database for MF6 (s26). BLAST's multi-database syntax sums the effective
database sizes, so E-values stay comparable across both; appending MF6 to the
1.9 GB FASTA and rebuilding would have been equivalent but wasteful.

THRESHOLDS ARE PRE-DECLARED and copied verbatim from s24 so they cannot drift
between the two searches:
  tier1  pident >= 60 and qcovs >= 70
  tier2  pident >= 60 and 30 <= qcovs < 70
  tier3  40 <= pident < 60 and qcovs >= 70
  tier4  anything else passing E <= 1e-5

Coverage is never dropped. Four of these queries begin with a signal peptide or
a lipobox (002734 ...AALALLAACGGG..., 003843, 004285), which will align to the
signal peptide of any secreted protein; 004284_1434_1538 is proline-rich
(PGGEREEPGRPYSPWPRRP). The query alignment span is exported for every hit so a
match carried entirely by such a region is visible rather than silently counted.
"""
import csv
import glob
import os
import subprocess
import sys
from collections import defaultdict

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
OUT = f"{W}/results"
SEARCH = f"{W}/search_toxins"
os.makedirs(SEARCH, exist_ok=True)
BIN = f"{W}/envs/bakta/bin"
DB = f"{W}/search_rhs/all_pyrodigal.faa {W}/pyrodigal_mf6/MF6.faa"

FMT = ("6 qseqid sseqid pident length nident mismatch gapopen qstart qend "
       "sstart send evalue bitscore qcovhsp qcovs slen qlen")
# FMT.split() is ["6", "qseqid", ...]; drop ONLY the "6" and zip against the
# whole row. Zipping COLS against f[1:] shifts every column silently.
COLS = FMT.split()[1:]

QUERIES = ["MF6_001079", "MF6_002734", "MF6_003184", "MF6_003684_2514_2867",
           "MF6_003684_full", "MF6_003686", "MF6_003843",
           "MF6_004284_1434_1538", "MF6_004284_full", "MF6_004285",
           "MF6_004947", "MF6_006318"]


def sh(cmd):
    print(f"$ {cmd[:200]}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                       env={**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}",
                            "TMPDIR": "/mnt/LargeStorageNoBackup/Moshea/tmp"})
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        sys.exit(f"FAILED rc={r.returncode}")
    return r


# ------------------------------------------------------- 1. contig size index
print("=== indexing contig lengths ===", flush=True)
contig_len, contig_rank = {}, {}
fnas = sorted(glob.glob(f"{W}/genomes/*.fna")) + [f"{W}/annot/MF6/MF6.fna"]
for fna in fnas:
    acc = "MF6" if fna.endswith("annot/MF6/MF6.fna") else os.path.basename(fna)[:-4]
    sizes, name, n = {}, None, 0
    with open(fna) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    sizes[name] = n
                name, n = line[1:].strip().split()[0], 0
            else:
                n += len(line.strip())
    if name:
        sizes[name] = n
    for r, (c, L) in enumerate(sorted(sizes.items(), key=lambda x: -x[1]), 1):
        contig_len[(acc, c)] = L
        contig_rank[(acc, c)] = r
ALL_ACC = sorted({("MF6" if f.endswith("annot/MF6/MF6.fna")
                   else os.path.basename(f)[:-4]) for f in fnas})
N_GEN = len(ALL_ACC)
print(f"  genomes: {N_GEN}   contigs: {len(contig_len)}", flush=True)
if N_GEN != 772:
    sys.exit(f"expected 772 genomes (771 + MF6), got {N_GEN}")

# ------------------------------------------------------- 2. blastp
tsv = f"{SEARCH}/blastp_hits.tsv"
if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
    sh(f'blastp -query {W}/queries_toxins12.faa -db "{DB}" -outfmt "{FMT}" '
       f'-evalue 1e-5 -max_target_seqs 1000000 -num_threads 48 '
       f'-comp_based_stats 2 -out {tsv}')
nrow = sum(1 for _ in open(tsv))
print(f"\nblastp rows: {nrow:,}", flush=True)

# ------------------------------------------------------- 3. metadata
meta, c3 = {}, {}
with open(f"{OUT}/replicon_census.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        meta[r["accession"]] = dict(organism=r["organism_name"],
                                    architecture=r["architecture"])
with open(f"{OUT}/c3_calls_all_genomes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c3[r["accession"]] = r["c3_present"]
meta["MF6"] = dict(organism="Burkholderia sola MF6", architecture="3_large")
c3["MF6"] = "True"


# ------------------------------------------------------- 4. tiers (from s24)
def tier(pid, qcov):
    if pid >= 60 and qcov >= 70:
        return 1
    if pid >= 60 and qcov >= 30:
        return 2
    if 40 <= pid < 60 and qcov >= 70:
        return 3
    return 4


hits = defaultdict(list)
bad = 0
with open(tsv) as fh:
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
        acc, contig, gidx, s, e, strand = parts
        pid, qcov = float(rec["pident"]), float(rec["qcovs"])
        hits[rec["qseqid"]].append(dict(
            acc=acc, contig=contig, gene=gidx, start=int(s), end=int(e),
            strand=strand, pident=pid, qcovs=qcov,
            qcovhsp=float(rec["qcovhsp"]), qstart=int(rec["qstart"]),
            qend=int(rec["qend"]), bits=float(rec["bitscore"]),
            evalue=float(rec["evalue"]), slen=int(rec["slen"]),
            qlen=int(rec["qlen"]), tier=tier(pid, qcov),
            contig_len=contig_len.get((acc, contig), 0),
            contig_rank=contig_rank.get((acc, contig), 0)))
if bad:
    sys.exit(f"ERROR: {bad} malformed rows - do not interpret these results")
print(f"parsed cleanly: {sum(len(v) for v in hits.values()):,} hits, 0 malformed",
      flush=True)

# ------------------------------------------------------- 5. tables
rows, summary = [], []
for q in QUERIES:
    hs = hits.get(q, [])
    by_acc = defaultdict(list)
    for h in hs:
        by_acc[h["acc"]].append(h)
    n_t1 = sum(1 for v in by_acc.values() if any(x["tier"] == 1 for x in v))
    n_t2 = sum(1 for v in by_acc.values()
               if not any(x["tier"] == 1 for x in v)
               and any(x["tier"] == 2 for x in v))
    n_t3 = sum(1 for v in by_acc.values()
               if not any(x["tier"] in (1, 2) for x in v)
               and any(x["tier"] == 3 for x in v))
    print(f"\n=== {q} ===", flush=True)
    print(f"  any hit E<=1e-5 : {len(by_acc)}/{N_GEN}", flush=True)
    print(f"  TIER1           : {n_t1}/{N_GEN}", flush=True)
    print(f"  TIER2 only      : {n_t2}", flush=True)
    print(f"  TIER3 only      : {n_t3}", flush=True)
    summary.append(dict(query=q, n_genomes_any_hit=len(by_acc), n_tier1=n_t1,
                        n_tier2_only=n_t2, n_tier3_only=n_t3,
                        n_genomes_total=N_GEN))
    for acc in ALL_ACC:
        v = by_acc.get(acc, [])
        best = max(v, key=lambda x: x["bits"]) if v else None
        rows.append(dict(
            query=q, accession=acc,
            organism=meta.get(acc, {}).get("organism", ""),
            architecture=meta.get(acc, {}).get("architecture", ""),
            pC3_present=c3.get(acc, ""),
            n_hits_any=len(v),
            n_hits_tier1=sum(1 for x in v if x["tier"] == 1),
            best_tier=best["tier"] if best else "",
            best_pident=f"{best['pident']:.1f}" if best else "",
            best_qcovs=f"{best['qcovs']:.1f}" if best else "",
            best_qstart=best["qstart"] if best else "",
            best_qend=best["qend"] if best else "",
            best_bitscore=f"{best['bits']:.0f}" if best else "",
            best_evalue=f"{best['evalue']:.1e}" if best else "",
            best_contig=best["contig"] if best else "",
            best_contig_len=best["contig_len"] if best else "",
            best_contig_rank=best["contig_rank"] if best else "",
            best_subject_len=best["slen"] if best else "",
            best_start=best["start"] if best else "",
            best_end=best["end"] if best else ""))

with open(f"{OUT}/toxin12_search_per_genome.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(f"{OUT}/toxin12_search_summary.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0]), delimiter="\t")
    w.writeheader(); w.writerows(summary)
with open(f"{OUT}/toxin12_search_all_hits.tsv", "w") as fh:
    cols = ["query", "acc", "contig", "gene", "start", "end", "strand", "pident",
            "qcovs", "qcovhsp", "qstart", "qend", "bits", "evalue", "slen",
            "qlen", "tier", "contig_len", "contig_rank"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for q in QUERIES:
        for h in sorted(hits.get(q, []), key=lambda x: -x["bits"]):
            w.writerow({**h, "query": q})

print("\nwrote toxin12_search_{per_genome,summary,all_hits}.tsv", flush=True)
