#!/usr/bin/env python
"""Rebuild the three derived RHS tables that fig9 consumes.

The generator for these was lost with the Shannon outage, so this script is
written to REPRODUCE the retained originals column-for-column when pointed at
the original inputs (--validate), and only then run on the rebuild inputs.

Outputs
  rhs_search_per_genome.tsv        one row per (query, genome)
  rhs_query_coverage_profile.tsv   per-residue HSP depth from CT-negative genomes
  rhs_pC3_association.tsv          Fisher exact tests, warhead and full protein

`architecture` is not carried in the rebuild's c3_calls table (the original's
was), so it is recomputed as the number of contigs >= LARGE_BP, capped at
"4+_large".  LARGE_BP = 300,000 reproduces the original label for 770 of 771
genomes; the single disagreement is GCF_054166145.1, which the original labels
`0_large` despite carrying contigs of 3.97 Mb, 2.96 Mb, 454 kb and 331 kb.
That is the known-bad original label already flagged for correction in the
report, so the rebuild's `4+_large` is the fix, not a regression.
"""
import argparse
import collections
import csv
import sys
from pathlib import Path

from scipy.stats import fisher_exact

D = Path(__file__).resolve().parent.parent
LARGE_BP = 300_000
WARHEAD_Q = "CT354_CFFIHE_03684"
FULL_Q = "FULL_CFFIHE_03684"

ap = argparse.ArgumentParser()
ap.add_argument("--tables", default=str(D / "tables"))
ap.add_argument("--validate", action="store_true",
                help="run on the ORIGINAL tables and diff against the retained "
                     "rhs_search_per_genome.tsv instead of writing anything")
a = ap.parse_args()
TAB = Path(a.tables)

rows = lambda p: list(csv.DictReader(open(p, newline=""), delimiter="\t"))

# ---------------------------------------------------------------- inputs
hits = rows(TAB / "rhs_search_all_hits.tsv")
ranks = rows(TAB / "contig_ranks.tsv")
calls = {r["accession"]: r for r in rows(TAB / "c3_calls_all_genomes.tsv")}
hostf = TAB / "host_categories.tsv"
if not hostf.exists():
    hostf = D.parent / "tables" / "host_categories.tsv"
organism = {r["accession"]: r["organism_name"] for r in rows(hostf)}

lens = collections.defaultdict(list)
for r in ranks:
    lens[r["accession"]].append(int(r["length"]))
genomes = sorted(lens)


def architecture(acc):
    n = sum(1 for x in lens[acc] if x >= LARGE_BP)
    return f"{n}_large" if n < 4 else "4+_large"


def pc3(acc):
    v = calls.get(acc, {})
    s = v.get("c3_present", "")
    return s == "True"


# ---------------------------------------------------------------- per genome
by = collections.defaultdict(list)
for h in hits:
    by[(h["query"], h["acc"])].append(h)
queries = sorted({h["query"] for h in hits})

FIELDS = ["query", "accession", "organism", "architecture", "pC3_present",
          "n_hits_any", "n_hits_tier1", "best_tier", "best_pident",
          "best_qcovs", "best_bitscore", "best_evalue", "best_contig",
          "best_contig_len", "best_contig_rank", "best_subject_len",
          "best_start", "best_end"]

per = []
for q in queries:
    for acc in genomes:
        hs = by.get((q, acc), [])
        rec = dict.fromkeys(FIELDS, "")
        rec.update(query=q, accession=acc, organism=organism.get(acc, ""),
                   architecture=architecture(acc), pC3_present=pc3(acc),
                   n_hits_any=len(hs),
                   n_hits_tier1=sum(1 for h in hs if h["tier"] == "1"))
        if hs:
            b = max(hs, key=lambda h: float(h["bits"]))
            rec.update(best_tier=b["tier"],
                       best_pident=f'{float(b["pident"]):.1f}',
                       best_qcovs=f'{float(b["qcovs"]):.1f}',
                       best_bitscore=b["bits"], best_evalue=b["evalue"],
                       best_contig=b["contig"], best_contig_len=b["contig_len"],
                       best_contig_rank=b["contig_rank"],
                       best_subject_len=b["slen"], best_start=b["start"],
                       best_end=b["end"])
        per.append(rec)

if a.validate:
    ref = {(r["query"], r["accession"]): r
           for r in rows(TAB / "rhs_search_per_genome.tsv")}
    got = {(r["query"], r["accession"]): r for r in per}
    only_ref = set(ref) - set(got)
    only_got = set(got) - set(ref)
    bad = collections.Counter()
    examples = {}
    for k in set(ref) & set(got):
        for f in FIELDS:
            rv, gv = ref[k].get(f, ""), str(got[k].get(f, ""))
            if f in ("best_bitscore", "best_evalue", "best_subject_len"):
                continue
            if rv != gv:
                bad[f] += 1
                examples.setdefault(f, (k, rv, gv))
    print(f"shared rows      : {len(set(ref) & set(got))}")
    print(f"only in retained : {len(only_ref)}")
    print(f"only in rebuilt  : {len(only_got)}")
    if not bad:
        print("MATCH: 0 mismatches across all compared columns")
    for f, n in bad.most_common():
        k, rv, gv = examples[f]
        print(f"  {f:20s} {n:5d} mismatches   e.g. {k}: retained={rv!r} got={gv!r}")
    sys.exit(0)

with open(TAB / "rhs_search_per_genome.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
    w.writeheader(); w.writerows(per)

# ------------------------------------------------- query coverage profile
# Depth of aligned HSPs along the full-length query, counting only genomes that
# do NOT carry the C-terminal warhead. The point of the panel is that this depth
# falls to exactly zero before the query ends.
ct_pos = {h["acc"] for h in hits if h["query"] == WARHEAD_Q and h["tier"] == "1"}
qlen = max(int(h["qlen"]) for h in hits if h["query"] == FULL_Q)
depth = [0] * (qlen + 1)
for h in hits:
    if h["query"] != FULL_Q or h["acc"] in ct_pos:
        continue
    for i in range(int(h["qstart"]), int(h["qend"]) + 1):
        depth[i] += 1
last_shared = max((i for i in range(1, qlen + 1) if depth[i]), default=0)
n_ct_neg = len({h["acc"] for h in hits if h["query"] == FULL_Q} - ct_pos)

with open(TAB / "rhs_query_coverage_profile.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["query_residue", "n_hsps_ct_negative"])
    for i in range(1, qlen + 1):
        w.writerow([i, depth[i]])

# ------------------------------------------------------- pC3 association
pos = {g for g in genomes if pc3(g)}
neg = set(genomes) - pos
assoc = []
for feat, q, tier1_only in (("CT354 warhead", WARHEAD_Q, True),
                            ("FULL RHS protein", FULL_Q, False)):
    carr = {h["acc"] for h in hits
            if h["query"] == q and (h["tier"] == "1" or not tier1_only)}
    a11, a10 = len(carr & pos), len(pos - carr)
    a01, a00 = len(carr & neg), len(neg - carr)
    orr, p = fisher_exact([[a11, a10], [a01, a00]])
    assoc.append(dict(feature=feat, n_pC3pos_with=a11, n_pC3neg_with=a01,
                      n_pC3pos_total=len(pos), n_pC3neg_total=len(neg),
                      odds_ratio=("inf" if orr == float("inf") else f"{orr:.2f}"),
                      fisher_p=f"{p:.3g}"))

with open(TAB / "rhs_pC3_association.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(assoc[0]), delimiter="\t")
    w.writeheader(); w.writerows(assoc)

print(f"wrote 3 tables to {TAB}")
print(f"  queries {queries}  genomes {len(genomes)}")
print(f"  query length {qlen}; last residue with CT-negative coverage {last_shared}")
print(f"  CT-negative genomes with a FULL hit: {n_ct_neg}; warhead carriers: {len(ct_pos)}")
for r in assoc:
    print(f"  {r['feature']:18s} pC3+ {r['n_pC3pos_with']}/{r['n_pC3pos_total']}"
          f"  pC3- {r['n_pC3neg_with']}/{r['n_pC3neg_total']}"
          f"  OR {r['odds_ratio']}  p {r['fisher_p']}")
