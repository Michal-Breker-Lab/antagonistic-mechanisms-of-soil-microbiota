#!/usr/bin/env python3
"""
s24 - search all 771 Burkholderia genomes for the MF6 RHS toxin.

Three queries, deliberately kept separate because they answer different questions:

  CT354_CFFIHE_03684  the C-terminal 354 aa toxin domain (what Moshe pasted)
                      -> "who carries THIS warhead?"
  FULL_CFFIHE_03684   the whole 2,867 aa RHS protein
                      -> "who carries an RHS protein of this family?"
  DEGQ_CFFIHE_02734   DegQ serine endoprotease
                      -> included because the accession was named; expected to be
                         near-universal, and that expectation is itself the control.

Comparing CT vs FULL is the point: RHS cores are conserved across genera while the
C-terminal cassette swaps. If a genome hits FULL but not CT, it has the delivery
scaffold with a DIFFERENT toxin.

THRESHOLDS ARE PRE-DECLARED HERE, before any output is inspected:

  tier1  pident >= 60 and qcovs >= 70   <- the answer to "60% ortholog"
  tier2  pident >= 60 and 30 <= qcovs < 70   partial: high local identity, short match
  tier3  40 <= pident < 60 and qcovs >= 70   diverged but full-length
  tier4  anything else passing E <= 1e-5     reported as counts only

The coverage floor in tier1 is not decoration. The first 180 aa of CT354 are 45%
Gly+Ala with AGGGAG repeats; without a length requirement an unrelated Gly-rich
membrane protein clears 60% identity over a 50-residue window.

-max_target_seqs is set very high rather than to a small N: the parameter keeps
the first N hits above threshold, NOT the best N (Shah et al. 2019,
doi:10.1093/bioinformatics/bty833). NCBI patched the worst of this in 2.8.1 and we
run 2.17.0, but a high ceiling makes the point moot.
"""
import csv
import glob
import os
import subprocess
import sys
from collections import defaultdict

W = "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3"
PROT = f"{W}/pyrodigal_prot"
OUT = f"{W}/results"
SEARCH = f"{W}/search_rhs"
os.makedirs(SEARCH, exist_ok=True)
BIN = f"{W}/envs/bakta/bin"

FMT = ("6 qseqid sseqid pident length nident mismatch gapopen qstart qend "
       "sstart send evalue bitscore qcovhsp qcovs slen qlen")
# FMT.split() is ["6", "qseqid", "sseqid", ...]; drop only the "6" so COLS lines
# up with the full row, qseqid included.
COLS = FMT.split()[1:]


def sh(cmd):
    print(f"$ {cmd[:160]}", flush=True)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                       env={**os.environ, "PATH": f"{BIN}:{os.environ['PATH']}",
                            "TMPDIR": f"{W}/tmp/tmpdir"})
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        sys.exit(f"FAILED rc={r.returncode}")
    return r


# ------------------------------------------------------------------ 1. contig sizes
print("=== indexing contig lengths ===", flush=True)
contig_len = {}          # (acc, contig) -> length
contig_rank = {}         # (acc, contig) -> rank by size, 1 = largest
for fna in sorted(glob.glob(f"{W}/genomes/*.fna")):
    acc = os.path.basename(fna)[:-4]
    sizes = {}
    name, n = None, 0
    with open(fna) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    sizes[name] = n
                name = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name:
        sizes[name] = n
    for r, (c, L) in enumerate(sorted(sizes.items(), key=lambda x: -x[1]), 1):
        contig_len[(acc, c)] = L
        contig_rank[(acc, c)] = r
print(f"  contigs indexed: {len(contig_len)}", flush=True)

# ------------------------------------------------------------------ 2. build DB
faas = sorted(glob.glob(f"{PROT}/*.faa"))
print(f"\n=== proteomes: {len(faas)} ===", flush=True)
if len(faas) < 771:
    sys.exit(f"expected 771 proteomes, found {len(faas)} - pyrodigal still running?")

allfaa = f"{SEARCH}/all_pyrodigal.faa"
if not os.path.exists(allfaa) or os.path.getsize(allfaa) == 0:
    with open(allfaa, "w") as out:
        for f in faas:
            with open(f) as fh:
                out.write(fh.read())
n_prot = sum(1 for line in open(allfaa) if line.startswith(">"))
print(f"  total predicted proteins: {n_prot:,}", flush=True)

if not os.path.exists(f"{allfaa}.pdb"):
    sh(f"makeblastdb -in {allfaa} -dbtype prot -out {allfaa} -title burk771")

# ------------------------------------------------------------------ 3. blastp
tsv = f"{SEARCH}/blastp_hits.tsv"
if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
    sh(f'blastp -query {W}/queries.faa -db {allfaa} -outfmt "{FMT}" '
       f'-evalue 1e-5 -max_target_seqs 1000000 -num_threads 48 '
       f'-comp_based_stats 2 -out {tsv}')
print(f"\nblastp rows: {sum(1 for _ in open(tsv)):,}", flush=True)

# ------------------------------------------------------------------ 4. metadata
meta = {}
with open(f"{OUT}/replicon_census.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        meta[r["accession"]] = dict(organism=r["organism_name"],
                                    architecture=r["architecture"])
c3 = {}
with open(f"{OUT}/c3_calls_all_genomes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c3[r["accession"]] = r["c3_present"]

ALL_ACC = sorted({os.path.basename(f)[:-4] for f in faas})

# ------------------------------------------------------------------ 5. tiers
def tier(pid, qcov):
    if pid >= 60 and qcov >= 70:
        return 1
    if pid >= 60 and qcov >= 30:
        return 2
    if 40 <= pid < 60 and qcov >= 70:
        return 3
    return 4


hits = defaultdict(list)
with open(tsv) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rec = dict(zip(COLS, f))
        q = rec["qseqid"]
        acc, contig, gidx, s, e, strand = rec["sseqid"].split("|")
        pid = float(rec["pident"]); qcov = float(rec["qcovs"])
        hits[q].append(dict(
            acc=acc, contig=contig, gene=gidx, start=int(s), end=int(e), strand=strand,
            pident=pid, qcovs=qcov, qcovhsp=float(rec["qcovhsp"]),
            bits=float(rec["bitscore"]), evalue=float(rec["evalue"]),
            slen=int(rec["slen"]), tier=tier(pid, qcov),
            contig_len=contig_len.get((acc, contig), 0),
            contig_rank=contig_rank.get((acc, contig), 0)))

rows = []
summary = []
for q in ["CT354_CFFIHE_03684", "FULL_CFFIHE_03684", "DEGQ_CFFIHE_02734"]:
    hs = hits.get(q, [])
    by_acc = defaultdict(list)
    for h in hs:
        by_acc[h["acc"]].append(h)
    n_t1 = sum(1 for a, v in by_acc.items() if any(x["tier"] == 1 for x in v))
    n_t2 = sum(1 for a, v in by_acc.items()
               if not any(x["tier"] == 1 for x in v) and any(x["tier"] == 2 for x in v))
    n_t3 = sum(1 for a, v in by_acc.items()
               if not any(x["tier"] in (1, 2) for x in v) and any(x["tier"] == 3 for x in v))
    n_any = len(by_acc)
    print(f"\n=== {q} ===", flush=True)
    print(f"  genomes with >=1 hit at E<=1e-5 : {n_any}/771", flush=True)
    print(f"  TIER1 >=60% id & >=70% cov      : {n_t1}/771", flush=True)
    print(f"  TIER2 >=60% id, 30-70% cov      : {n_t2}", flush=True)
    print(f"  TIER3 40-60% id, >=70% cov      : {n_t3}", flush=True)
    summary.append(dict(query=q, n_genomes_any_hit=n_any, n_tier1=n_t1,
                        n_tier2_only=n_t2, n_tier3_only=n_t3, n_genomes_total=771))
    for acc in ALL_ACC:
        v = by_acc.get(acc, [])
        best = max(v, key=lambda x: x["bits"]) if v else None
        t1 = [x for x in v if x["tier"] == 1]
        rows.append(dict(
            query=q, accession=acc,
            organism=meta.get(acc, {}).get("organism", ""),
            architecture=meta.get(acc, {}).get("architecture", ""),
            pC3_present=c3.get(acc, ""),
            n_hits_any=len(v), n_hits_tier1=len(t1),
            best_tier=best["tier"] if best else "",
            best_pident=f"{best['pident']:.1f}" if best else "",
            best_qcovs=f"{best['qcovs']:.1f}" if best else "",
            best_bitscore=f"{best['bits']:.0f}" if best else "",
            best_evalue=f"{best['evalue']:.1e}" if best else "",
            best_contig=best["contig"] if best else "",
            best_contig_len=best["contig_len"] if best else "",
            best_contig_rank=best["contig_rank"] if best else "",
            best_subject_len=best["slen"] if best else "",
            best_start=best["start"] if best else "",
            best_end=best["end"] if best else ""))

with open(f"{OUT}/rhs_search_per_genome.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(f"{OUT}/rhs_search_summary.tsv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(summary)

# every individual hit, unfiltered, for auditing
with open(f"{OUT}/rhs_search_all_hits.tsv", "w") as fh:
    cols = ["query", "acc", "contig", "gene", "start", "end", "strand", "pident",
            "qcovs", "qcovhsp", "bits", "evalue", "slen", "tier",
            "contig_len", "contig_rank"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for q, hs in hits.items():
        for h in sorted(hs, key=lambda x: -x["bits"]):
            w.writerow({**h, "query": q})

print("\nwrote rhs_search_per_genome.tsv, rhs_search_summary.tsv, "
      "rhs_search_all_hits.tsv", flush=True)
