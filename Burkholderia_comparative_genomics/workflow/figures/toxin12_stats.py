#!/usr/bin/env python3
"""Statistics for the 12-query / 10-locus toxin search across 773 genomes.

The 12 sequences Moshe supplied are 10 LOCI. Two are searched at two scales
because that contrast is the point for polymorphic toxins:

    MF6_003684  RHS repeat protein, 2867 aa  +  its C-terminal 354 aa
    MF6_004284  Rhs family protein, 1538 aa  +  its C-terminal 104 aa

Counting them as 12 would double-count those two loci in every test, so the
locus-level statistics use the FULL-LENGTH query and the C-terminal query is
reported as a within-locus contrast.

Locus tags are PGAP (Andrei's annotation) because that is the numbering the
queries came from; the Bakta CFFIHE_ tag is carried as a cross-reference
because every other table in this project is keyed on it.

One locus is not searched by blastp. MF6_003686 is a 120-aa ORF that Bakta
missed entirely and that pyrodigal truncates to 78 aa by choosing an internal
start codon -- MF6 scores only 65% coverage against its OWN gene, below the
tier-1 floor. A protein-level search therefore cannot answer this locus at all,
so it is searched with tblastn against the assemblies instead. That is a
DIFFERENT UNIT OF SEARCH (six-frame translation of DNA, not predicted proteins)
and is flagged as such in the output; qcovhsp is used in place of qcovs.

Tests, all pre-declared:
  * Fisher exact per locus for association with pC3 presence, over all 773
    genomes only (MF6 is the query source and would bias every cell), then
    Benjamini-Hochberg across the 10 loci.
  * Paralogy diagnostic: a locus where most genomes carry several hits is a gene
    FAMILY, not an ortholog, and its map must not be read as effector presence.
  * Composition-bias audit: tier-1 hits whose aligned query span falls entirely
    within the first 40 residues are carried by a signal peptide / lipobox, not
    by the body of the protein.
"""
import csv
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

D = Path(__file__).resolve().parent.parent
TAB = D / "tables"

# ---------------------------------------------------------------- locus table
# PGAP tag, Bakta tag, replicon, length, PGAP product, eggNOG description
LOCI = [
    ("MF6_001079", "CFFIHE_01081", "chr1", 229, "hypothetical protein",
     "(no eggNOG annotation)"),
    ("MF6_002734", "CFFIHE_02739", "chr1", 379, "SGNH/GDSL hydrolase family protein",
     "GDSL-like Lipase/Acylhydrolase (Lipase_GDSL)"),
    ("MF6_003184", "CFFIHE_03187", "chr1", 854, "T6SS effector BTH_I2691 family protein",
     "Band 7 protein"),
    ("MF6_003684", "CFFIHE_03684", "pC3", 2867, "RHS repeat domain-containing protein",
     "self proteolysis"),
    ("MF6_003686", "(not in Bakta)", "pC3", 120, "hypothetical protein",
     "(no eggNOG annotation)"),
    ("MF6_003843", "CFFIHE_03842", "pC3", 565, "M4 family metallopeptidase",
     "Thermolysin metallopeptidase (Peptidase_M4), ko:K20273"),
    ("MF6_004284", "CFFIHE_04281", "pC3", 1538, "RHS repeat-associated core domain protein",
     "COG3209 Rhs family protein (RHS_repeat)"),
    ("MF6_004285", "CFFIHE_04283", "pC3", 326, "hypothetical protein",
     "(no eggNOG annotation)"),
    ("MF6_004947", "CFFIHE_04987", "chr2", 570, "M36 family metallopeptidase",
     "Fungalysin metallopeptidase M36, ko:K20274"),
    ("MF6_006318", "CFFIHE_06334", "chr2", 364, "triacylglycerol lipase",
     "Triacylglycerol lipase (Abhydrolase_1), ko:K01046"),
]
# locus -> (full-length query id, C-terminal query id or None)
QMAP = {
    "MF6_001079": ("MF6_001079", None),
    "MF6_002734": ("MF6_002734", None),
    "MF6_003184": ("MF6_003184", None),
    "MF6_003684": ("MF6_003684_full", "MF6_003684_2514_2867"),
    "MF6_003686": ("MF6_003686", None),          # tblastn, see module docstring
    "MF6_003843": ("MF6_003843", None),
    "MF6_004284": ("MF6_004284_full", "MF6_004284_1434_1538"),
    "MF6_004285": ("MF6_004285", None),
    "MF6_004947": ("MF6_004947", None),
    "MF6_006318": ("MF6_006318", None),
}
TBLASTN_LOCI = {"MF6_003686"}


def tier(pid, qcov):
    """Copied verbatim from s24_search.py -- must not drift between searches."""
    if pid >= 60 and qcov >= 70:
        return 1
    if pid >= 60 and qcov >= 30:
        return 2
    if 40 <= pid < 60 and qcov >= 70:
        return 3
    return 4


# ---------------------------------------------------------------- blastp table
pg = defaultdict(dict)          # query -> acc -> row
c3, organism = {}, {}
with open(TAB / "toxin12_search_per_genome.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        pg[r["query"]][r["accession"]] = r
        c3[r["accession"]] = r["pC3_present"]
        organism[r["accession"]] = r["organism"]
ACCS = sorted(c3)
# D14: MF6 and MF7 are INCLUDED in the association tests. The original excluded
# MF6 as the query source (every query is an MF6 locus, so MF6 self-hits at
# 100%/100%); the rebuild counts them, with the circularity stated in Methods
# rather than hidden by an undocumented exclusion. TESTSET is therefore all
# genomes, and n_carriers_773 replaces the original's n_carriers_771.
TESTSET = ACCS
_q = [a for a in ACCS if a in ("MF6", "MF7")]
print(f"genomes: {len(ACCS)} (association tested over all {len(TESTSET)}; "
      f"query-source genomes included: {', '.join(_q) or 'none'})")

# ---------------------------------------------------------------- tblastn table
tb_cols = ("qseqid sseqid pident length qstart qend sstart send evalue bitscore "
           "qcovhsp qlen").split()
tb = defaultdict(list)
with open(TAB / "toxin12_tblastn_003686.tsv") as fh:
    for line in fh:
        r = dict(zip(tb_cols, line.rstrip("\n").split("\t")))
        acc = r["sseqid"].split("|", 1)[0]
        tb[acc].append(r)

# ---------------------------------------------------------------- presence calls
# carriers[locus][acc] = True if tier 1
carriers, best_of = {}, {}
for tag, _, _, _, _, _ in LOCI:
    full_q = QMAP[tag][0]
    carriers[tag], best_of[tag] = {}, {}
    for acc in ACCS:
        if tag in TBLASTN_LOCI:
            hs = tb.get(acc, [])
            if hs:
                b = max(hs, key=lambda h: float(h["bitscore"]))
                pid, qcov = float(b["pident"]), float(b["qcovhsp"])
                best_of[tag][acc] = (pid, qcov, len(hs))
                carriers[tag][acc] = tier(pid, qcov) == 1
            else:
                best_of[tag][acc] = (0.0, 0.0, 0)
                carriers[tag][acc] = False
        else:
            r = pg[full_q][acc]
            n = int(r["n_hits_any"])
            if n:
                best_of[tag][acc] = (float(r["best_pident"]),
                                     float(r["best_qcovs"]), n)
            else:
                best_of[tag][acc] = (0.0, 0.0, 0)
            carriers[tag][acc] = int(r["n_hits_tier1"]) > 0

# C-terminal queries, reported as a within-locus contrast (never as a locus)
ct_carriers = {}
for tag, (full_q, ct_q) in QMAP.items():
    if ct_q:
        ct_carriers[tag] = {a: int(pg[ct_q][a]["n_hits_tier1"]) > 0 for a in ACCS}

# ---------------------------------------------------------------- Fisher + BH
rows, pvals = [], []
for tag, bakta, rep, ln, prod, egg in LOCI:
    a = sum(1 for x in TESTSET if carriers[tag][x] and c3[x] == "True")     # carrier, pC3+
    b = sum(1 for x in TESTSET if carriers[tag][x] and c3[x] != "True")     # carrier, pC3-
    cc = sum(1 for x in TESTSET if not carriers[tag][x] and c3[x] == "True")
    d = sum(1 for x in TESTSET if not carriers[tag][x] and c3[x] != "True")
    odds, p = fisher_exact([[a, b], [cc, d]])
    pvals.append(p)
    rows.append(dict(locus=tag, bakta_tag=bakta, replicon=rep, length_aa=ln,
                     pgap_product=prod, eggnog=egg,
                     search="tblastn" if tag in TBLASTN_LOCI else "blastp",
                     n_carriers_773=a + b, n_carrier_pC3pos=a, n_carrier_pC3neg=b,
                     n_noncarrier_pC3pos=cc, n_noncarrier_pC3neg=d,
                     odds_ratio=f"{odds:.3g}", p_raw=f"{p:.3g}"))
rej, q, _, _ = multipletests(pvals, method="fdr_bh")
for r, qq, rj in zip(rows, q, rej):
    r["p_bh"] = f"{qq:.3g}"
    r["significant_bh_0.05"] = bool(rj)

# ---------------------------------------------------------------- paralogy
for r in rows:
    tag = r["locus"]
    ns = sorted(best_of[tag][a][2] for a in TESTSET)
    hit_g = [n for n in ns if n]
    r["median_hits_in_hit_genomes"] = (
        sorted(hit_g)[len(hit_g) // 2] if hit_g else 0)
    r["max_hits_in_a_genome"] = max(ns) if ns else 0
    r["frac_hit_genomes_multi"] = (
        f"{sum(1 for n in hit_g if n >= 2) / len(hit_g):.2f}" if hit_g else "NA")
    # a locus is low-information if it is both near-universal and multi-copy
    r["low_information"] = bool(
        r["n_carriers_773"] > 0.5 * len(TESTSET)
        and hit_g and sum(1 for n in hit_g if n >= 2) / len(hit_g) > 0.20)

# ---------------------------------------------------------------- comp-bias audit
span_bad = defaultdict(int)
with open(TAB / "toxin12_search_all_hits.tsv") as fh:
    for h in csv.DictReader(fh, delimiter="\t"):
        if int(h["tier"]) == 1 and int(h["qend"]) <= 40:
            span_bad[h["query"]] += 1
for r in rows:
    fq = QMAP[r["locus"]][0]
    r["tier1_hits_in_first40aa_only"] = span_bad.get(fq, 0)

with open(TAB / "toxin12_locus_summary.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ---------------------------------------------------------------- report
print("\n" + "=" * 108)
print(f"{'locus':<12}{'rep':<6}{'aa':>6}{'carriers/773':>14}{'pC3+':>6}{'pC3-':>6}"
      f"{'OR':>9}{'p_BH':>10}  {'multi':>6} {'lowinfo':>8}  product")
print("=" * 108)
for r in rows:
    print(f"{r['locus']:<12}{r['replicon']:<6}{r['length_aa']:>6}"
          f"{r['n_carriers_773']:>14}{r['n_carrier_pC3pos']:>6}"
          f"{r['n_carrier_pC3neg']:>6}{r['odds_ratio']:>9}{r['p_bh']:>10}"
          f"  {r['frac_hit_genomes_multi']:>6} {str(r['low_information']):>8}"
          f"  {r['pgap_product'][:34]}")

print(f"\nwithin-locus contrast: full-length vs C-terminal tip (tier-1, of {len(ACCS)})")
for tag, (fq, ct) in QMAP.items():
    if ct:
        nf = sum(1 for a in ACCS if carriers[tag][a])
        nc = sum(1 for a in ACCS if ct_carriers[tag][a])
        nany = sum(1 for a in ACCS if best_of[tag][a][2] > 0)
        print(f"  {tag}: any hit {nany:>3}   full-length tier1 {nf:>3}   "
              f"C-terminal tier1 {nc:>3}")

print("\ncomposition-bias audit (tier-1 hits carried only by residues 1-40):")
for r in rows:
    print(f"  {r['locus']}: {r['tier1_hits_in_first40aa_only']}")

# carrier sets, for the figure and the phylogenetic test
with open(TAB / "toxin12_carriers.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["locus", "accession", "organism", "pC3_present", "tier1_carrier",
                "best_pident", "best_qcov", "n_hits"])
    for tag, _, _, _, _, _ in LOCI:
        for acc in ACCS:
            pid, qc, n = best_of[tag][acc]
            w.writerow([tag, acc, organism.get(acc, ""), c3[acc],
                        carriers[tag][acc], f"{pid:.1f}", f"{qc:.1f}", n])
    for tag, cm in ct_carriers.items():
        for acc in ACCS:
            r = pg[QMAP[tag][1]][acc]
            w.writerow([tag + "_Cterm", acc, organism.get(acc, ""), c3[acc],
                        cm[acc], r["best_pident"] or "0.0",
                        r["best_qcovs"] or "0.0", r["n_hits_any"]])

print("\nwrote toxin12_locus_summary.tsv, toxin12_carriers.tsv")
