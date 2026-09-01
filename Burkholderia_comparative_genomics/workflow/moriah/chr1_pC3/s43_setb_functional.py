#!/usr/bin/env python3
"""Functional characterisation of the Set B clade cores -> the fig7 tables.

Rebuilds the five tables the original (lost) generator produced, from:
  - the family tables written by s42 (mmseqs easy-cluster --min-seq-id 0.8
    -c 0.8 --cov-mode 1; those settings were recovered by sweeping until the
    retained pC3 counts reproduced exactly -- 1,382 families / 833 core / 549
    accessory);
  - Bakta's DbXrefs column, which already carries COG id, COG category, KEGG KO
    and product for every CDS, so nothing is re-annotated;
  - the genus-wide pC3 pangenome's soft-core partition, bridged through MF6
    locus tags, to split the clade core into genus-derived and clade-specific.

Core is SINGLE-COPY core: present in every genome with exactly one gene each.
The 10 pC3 families that are present in all genomes but multi-copy in at least
one are counted as accessory, which is what makes the totals 833 + 549 = 1,382
rather than 843 + 539.

Family-level annotation is by majority vote across members (report section 6.4); counting
per gene would multiply every count by the genome count.

Statistics are descriptive. The two sets being contrasted are two replicons of
the same genomes, not independent samples, so Fisher's null of random
gene-to-replicon assignment is not a hypothesis anyone holds; log2 odds ratios
with Haldane-Anscombe-corrected 95% CIs lead, and BH q values are reported only
to mark which contrasts are robust.
"""
import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import fisher_exact
import numpy as np

COG_NAMES = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control and mitosis",
    "E": "Amino acid metabolism and transport", "F": "Nucleotide metabolism and transport",
    "G": "Carbohydrate metabolism and transport", "H": "Coenzyme metabolism",
    "I": "Lipid metabolism", "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown",
    "T": "Signal transduction mechanisms", "U": "Intracellular trafficking and secretion",
    "V": "Defence mechanisms", "W": "Extracellular structures",
    "X": "Mobilome: prophages, transposons", "Y": "Nuclear structure",
    "Z": "Cytoskeleton", "-": "no COG assigned",
}


def read_tsv(p):
    return list(csv.DictReader(open(p, newline=""), delimiter="\t"))


def bakta_annotation(annot_dir, accessions):
    """locus_tag -> dict(cog_id, cog_cats, kos, product), from Bakta's TSV."""
    ann = {}
    for acc in accessions:
        p = Path(annot_dir) / acc / f"{acc}.tsv"
        for line in open(p):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[1] != "cds":
                continue
            xr = [x.strip() for x in f[8].split(",")]
            # Bakta writes both "COG:COG#####" (the orthologue) and "COG:<letters>"
            # (its functional category) into the same field.
            cog_id = next((x[4:] for x in xr if x.startswith("COG:COG")), "")
            cats = next((x[4:] for x in xr
                         if x.startswith("COG:") and not x.startswith("COG:COG")), "")
            kos = ";".join(sorted({x[6:] for x in xr if x.startswith("KEGG:")}))
            ann[f[5]] = dict(cog_id=cog_id, cog_cats=cats, kos=kos, product=f[7])
    return ann


def majority(values, n_members):
    """Value held by a STRICT majority of the family's members, else "".

    Not "most common non-empty": with 5 members, a COG carried by one member and
    absent from four is not the family's annotation. Requiring > n/2 is what
    reproduces the retained coverage figures (239/833 = 28.7% of the pC3 clade
    core carrying a COG category, not 30.9%).
    """
    c = Counter(v for v in values if v)
    if not c:
        return ""
    val, n = c.most_common(1)[0]
    ties = sorted(k for k, m in c.items() if m == n)
    return ties[0] if n * 2 > n_members else ""


def bucket(product):
    p = product.lower()
    if "hypothetical" in p:
        return "hypothetical protein"
    if re.search(r"\bduf\d+", p) or "uncharacteri" in p:
        return "DUF / uncharacterised domain"
    if "domain-containing protein" in p:
        return "domain-containing, no specific function"
    return "named product"


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        prev = min(prev, p[i] * n / (n - rank + 1))
        q[i] = prev
    return q


def enrich(cats_a, n_a, cats_b, n_b, set_a, set_b, contrast, mode):
    rows = []
    for cat in sorted(set(cats_a) | set(cats_b) | {"-"}):
        a, b = cats_a.get(cat, 0), cats_b.get(cat, 0)
        # Haldane-Anscombe: +0.5 to every cell, so a zero cell gives a finite OR
        t = [[a + 0.5, n_a - a + 0.5], [b + 0.5, n_b - b + 0.5]]
        lor = np.log2((t[0][0] * t[1][1]) / (t[0][1] * t[1][0]))
        se = np.sqrt(sum(1 / x for x in (t[0][0], t[0][1], t[1][0], t[1][1]))) / np.log(2)
        _, p = fisher_exact([[a, n_a - a], [b, n_b - b]])
        rows.append(dict(mode=mode, contrast=contrast, cog_category=cat,
                         cog_name=COG_NAMES.get(cat, cat),
                         set_a=set_a, n_a=a, total_a=n_a, pct_a=round(100 * a / n_a, 2),
                         set_b=set_b, n_b=b, total_b=n_b, pct_b=round(100 * b / n_b, 2),
                         log2_OR=round(lor, 3), ci_low=round(lor - 1.96 * se, 3),
                         ci_high=round(lor + 1.96 * se, 3), p_fisher=p))
    for r, q in zip(rows, bh([r["p_fisher"] for r in rows])):
        r["q_BH"] = q
    return rows


ap = argparse.ArgumentParser()
ap.add_argument("--fam-pc3", required=True, type=Path)
ap.add_argument("--fam-chr1", required=True, type=Path)
ap.add_argument("--annot", required=True, type=Path)
ap.add_argument("--genus-panaroo-csv", type=Path,
                help="Panaroo gene_presence_absence.csv for the genus-wide pC3\n                     pangenome; preferred over the PPanGGOLiN partition")
ap.add_argument("--genus-softcore-frac", type=float, default=0.95)
ap.add_argument("--mf6-accession", default="MF6")
ap.add_argument("--genus-softcore", type=Path,
                help="PPanGGOLiN partitions/soft_core.txt for the genus-wide pC3 pangenome")
ap.add_argument("--genus-families", type=Path,
                help="PPanGGOLiN gene_families.tsv (col1 = family, col3 = Bakta locus tag)")
ap.add_argument("--outdir", required=True, type=Path)
ap.add_argument("--mf6-prefix", default="CFFIHE",
                help="MF6's Bakta locus-tag prefix")
ap.add_argument("--label", required=True, help="tag for this run, e.g. 5g or 6g")
a = ap.parse_args()
a.outdir.mkdir(parents=True, exist_ok=True)

fam_pc3 = read_tsv(a.fam_pc3)
fam_chr1 = read_tsv(a.fam_chr1)
accs = sorted({m.split("|")[0] for r in fam_pc3 for m in r["members"].split(";")})
ann = bakta_annotation(a.annot, accs)
print(f"{len(accs)} genomes, {len(ann):,} annotated CDS")

# ---- genus soft-core families, mapped to the Bakta locus tags they contain ---
# Bridge through MF6's locus tags ONLY (report section 6.4). A genus soft-core family
# spans ~254 genomes, so matching ANY member's tag would call almost every clade
# family genus-derived; the question is specifically which of MF6's pC3 genes sit
# in a genus-wide core family.
mf6_prefix = a.mf6_prefix + "_"
genus_tags = set()
if a.genus_panaroo_csv:
    # Panaroo's gene_presence_absence.csv holds one column per genome carrying
    # that genome's locus tag(s) for the family, so MF6's tags come straight out
    # of the MF6 column -- no separate family->tag map needed.
    import io
    with open(a.genus_panaroo_csv, newline="") as fh:
        rd = csv.reader(fh)
        hdr = next(rd)
        gcols = list(range(3, len(hdr)))
        try:
            mf6col = hdr.index(a.mf6_accession)
        except ValueError:
            raise SystemExit(f"FAIL: no column {a.mf6_accession} in {a.genus_panaroo_csv}")
        n_soft = 0
        for row in rd:
            present = sum(1 for i in gcols if i < len(row) and row[i].strip())
            if present < a.genus_softcore_frac * len(gcols):
                continue
            n_soft += 1
            for t in row[mf6col].split(";"):
                t = t.strip()
                if t.startswith(mf6_prefix):
                    genus_tags.add(t)
    print(f"genus pC3 soft-core (Panaroo, >={a.genus_softcore_frac:.0%} of "
          f"{len(gcols)} genomes): {n_soft} families -> {len(genus_tags):,} MF6 locus tags")
else:
    softcore = {l.strip() for l in open(a.genus_softcore) if l.strip()}
    for line in open(a.genus_families):
        f = line.rstrip("\n").split("\t")
        if f[0] in softcore and len(f) > 2 and f[2].startswith(mf6_prefix):
            genus_tags.add(f[2])
    print(f"genus pC3 soft-core (PPanGGOLiN): {len(softcore)} families -> "
          f"{len(genus_tags):,} MF6 locus tags")


def annotate(rows, gene_set_of):
    out = []
    for r in rows:
        mem = [m.split("|", 1) for m in r["members"].split(";")]
        tags = [t for _, t in mem]
        recs = [ann[t] for t in tags if t in ann]
        gs = gene_set_of(r, tags)
        if gs is None:
            continue
        nm = len(mem)
        out.append(dict(
            gene_set=gs, family=r["family"],
            cog_id=majority((x["cog_id"] for x in recs), nm),
            cog_categories=majority((x["cog_cats"] for x in recs), nm),
            kos=majority((x["kos"] for x in recs), nm),
            n_members=len(mem), n_genomes=int(r["n_genomes"]),
            genomes=r["genomes"],
            # Product text is free-form and varies in wording between genomes, so a
            # strict majority rarely carries; take the most common instead. COG/KO
            # are controlled vocabularies and do take the strict rule.
            product=(Counter(x["product"] for x in recs if x["product"]).most_common(1)
                     or [("hypothetical protein", 0)])[0][0]))
    return out


def pc3_set(r, tags):
    if r["is_single_copy_core"] == "1":
        sub = "genus_derived" if any(t in genus_tags for t in tags) else "clade_specific"
        return f"pC3_clade_core/{sub}"
    return "pC3_clade_accessory"


def chr1_set(r, tags):
    return "chr1_clade_core" if r["is_single_copy_core"] == "1" else None


fams = annotate(fam_pc3, pc3_set) + annotate(fam_chr1, chr1_set)
with open(a.outdir / f"setB_functional_families_{a.label}.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["gene_set", "family", "cog_id", "cog_categories", "kos",
                            "n_members", "n_genomes", "genomes", "product"],
                       delimiter="\t")
    w.writeheader()
    w.writerows(fams)

SETS = ["pC3_clade_core", "pC3_clade_accessory", "chr1_clade_core",
        "pC3_core_genus_derived", "pC3_core_clade_specific"]


def members_of(name):
    if name == "pC3_core_genus_derived":
        return [f for f in fams if f["gene_set"] == "pC3_clade_core/genus_derived"]
    if name == "pC3_core_clade_specific":
        return [f for f in fams if f["gene_set"] == "pC3_clade_core/clade_specific"]
    return [f for f in fams if f["gene_set"].split("/")[0] == name]


# ---- COG profile (a family with several categories counts in each) ----------
prof = []
for s in SETS[:3]:
    ms = members_of(s)
    n = len(ms)
    c = Counter()
    for f in ms:
        cc = f["cog_categories"] or "-"
        for ch in (cc if cc != "-" else "-"):
            c[ch] += 1
    for cat in sorted(c):
        prof.append(dict(gene_set=s, cog_category=cat, cog_name=COG_NAMES.get(cat, cat),
                         n_families=c[cat], total=n, pct=round(100 * c[cat] / n, 3)))
with open(a.outdir / f"setB_cog_profile_{a.label}.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["gene_set", "cog_category", "cog_name", "n_families",
                            "total", "pct"], delimiter="\t")
    w.writeheader(); w.writerows(prof)


def cat_counts(ms, annotated_only):
    ms = [f for f in ms if f["cog_categories"]] if annotated_only else ms
    c = Counter()
    for f in ms:
        cc = f["cog_categories"] or "-"
        for ch in (cc if cc != "-" else "-"):
            c[ch] += 1
    return c, len(ms)


enr = []
CONTRASTS = [("A_pC3core_vs_chr1core", "pC3_clade_core", "chr1_clade_core"),
             ("B_pC3core_vs_pC3acc", "pC3_clade_core", "pC3_clade_accessory"),
             ("C_genusderived_vs_cladespecific", "pC3_core_genus_derived",
              "pC3_core_clade_specific")]
for mode, aonly in (("all_families", False), ("COG_annotated_only", True)):
    for name, sa, sb in CONTRASTS:
        ca, na = cat_counts(members_of(sa), aonly)
        cb, nb = cat_counts(members_of(sb), aonly)
        enr += enrich(ca, na, cb, nb, sa, sb, name, mode)
with open(a.outdir / f"setB_cog_enrichment_{a.label}.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, ["mode", "contrast", "cog_category", "cog_name", "set_a",
                            "n_a", "total_a", "pct_a", "set_b", "n_b", "total_b",
                            "pct_b", "log2_OR", "ci_low", "ci_high", "p_fisher",
                            "q_BH"], delimiter="\t")
    w.writeheader(); w.writerows(enr)

# ---- annotation coverage and product buckets -------------------------------
cov, buck = [], []
for s in SETS:
    ms = members_of(s)
    n = len(ms)
    nc = sum(1 for f in ms if f["cog_categories"])
    nk = sum(1 for f in ms if f["kos"])
    nh = sum(1 for f in ms if "hypothetical" in f["product"].lower())
    cov.append(dict(gene_set=s, n_families=n, n_with_COG=nc,
                    pct_with_COG=round(100 * nc / n, 1), n_with_KO=nk,
                    pct_with_KO=round(100 * nk / n, 1), n_hypothetical=nh,
                    pct_hypothetical=round(100 * nh / n, 1)))
    if s in SETS[:3]:
        bc = Counter(bucket(f["product"]) for f in ms)
        for b in sorted(bc, key=lambda x: -bc[x]):
            buck.append(dict(gene_set=s, product_bucket=b, n_families=bc[b],
                             total=n, pct=round(100 * bc[b] / n, 2)))
for name, rows, cols in (
        (f"setB_annotation_coverage_{a.label}.tsv", cov,
         ["gene_set", "n_families", "n_with_COG", "pct_with_COG", "n_with_KO",
          "pct_with_KO", "n_hypothetical", "pct_hypothetical"]),
        (f"setB_product_buckets_{a.label}.tsv", buck,
         ["gene_set", "product_bucket", "n_families", "total", "pct"])):
    with open(a.outdir / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader(); w.writerows(rows)

summary = dict(setB_genomes=accs,
               **{s: len(members_of(s)) for s in SETS},
               genus_softcore_mf6_tags=len(genus_tags),
               coverage=cov)
(a.outdir / f"setB_functional_summary_{a.label}.json").write_text(
    json.dumps(summary, indent=2))

print(f"\n{'gene set':<28}{'families':>10}{'%COG':>8}{'%KO':>8}{'%hypo':>8}")
for c in cov:
    print(f"{c['gene_set']:<28}{c['n_families']:>10,}{c['pct_with_COG']:>8}"
          f"{c['pct_with_KO']:>8}{c['pct_hypothetical']:>8}")
print(f"\nwrote 5 tables + summary to {a.outdir} (label {a.label})")
