#!/usr/bin/env python3
"""Ward-cluster the day-2 core proteome and test each cluster for COG and KEGG
over-representation.
"""
from types import SimpleNamespace
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[2]

STRAIN_ORDER = ["MF6", "27D6", "34F7"]

COG_CATEGORIES = {
    "J": ("Translation, ribosomal structure and biogenesis", "Information storage"),
    "A": ("RNA processing and modification", "Information storage"),
    "K": ("Transcription", "Information storage"),
    "L": ("Replication, recombination and repair", "Information storage"),
    "B": ("Chromatin structure and dynamics", "Information storage"),
    "D": ("Cell cycle control, cell division, chromosome partitioning", "Cellular processes"),
    "Y": ("Nuclear structure", "Cellular processes"),
    "V": ("Defense mechanisms", "Cellular processes"),
    "T": ("Signal transduction mechanisms", "Cellular processes"),
    "M": ("Cell wall/membrane/envelope biogenesis", "Cellular processes"),
    "N": ("Cell motility", "Cellular processes"),
    "Z": ("Cytoskeleton", "Cellular processes"),
    "W": ("Extracellular structures", "Cellular processes"),
    "U": ("Intracellular trafficking, secretion, vesicular transport", "Cellular processes"),
    "O": ("Post-translational modification, protein turnover, chaperones", "Cellular processes"),
    "C": ("Energy production and conversion", "Metabolism"),
    "G": ("Carbohydrate transport and metabolism", "Metabolism"),
    "E": ("Amino acid transport and metabolism", "Metabolism"),
    "F": ("Nucleotide transport and metabolism", "Metabolism"),
    "H": ("Coenzyme transport and metabolism", "Metabolism"),
    "I": ("Lipid transport and metabolism", "Metabolism"),
    "P": ("Inorganic ion transport and metabolism", "Metabolism"),
    "Q": ("Secondary metabolites biosynthesis, transport, catabolism", "Metabolism"),
    "R": ("General function prediction only", "Poorly characterized"),
    "S": ("Function unknown", "Poorly characterized"),
    "X": ("Mobilome: prophages, transposons", "Poorly characterized"),
}

def kegg_known_set(path):
    """Every map id KEGG currently publishes, from the BRITE class table.

    An id in the annotation that is NOT here has been retired: eggNOG's KEGG
    assignment lags the current release, so ids KEGG has withdrawn keep arriving.
    ko01130 "Biosynthesis of antibiotics" is the case that motivated this - a
    retired global map with 196 proteins that no class filter could catch, because
    it has no class to filter on, and no name in kegg_pathway_names.tsv either, so
    it plotted as a bare id.  A term with no name and no definition cannot be
    interpreted, and it is the same superset it always was.
    """
    with open(path) as fh:
        return {r["pathway_id"] for r in
                csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")}


def kegg_drop_set(path, drop_non_bacterial, drop_overview):
    """ko##### ids to exclude, from resources/annotation/kegg_pathway_classes.tsv.

    non-bacterial : Organismal Systems / Human Diseases / Drug Development.  eggNOG
        reaches these through KOs shared between a bacterial enzyme and a
        eukaryote-specific pathway - that is how "PPAR signaling pathway" turns up
        in a Burkholderia.
    overview      : the "Global and overview maps" subclass.  ko01100 alone covers
        ~46% of the KEGG-annotated proteome; these are supersets, not pathways.

    Dropped BEFORE testing rather than flagged afterwards: left in, they enlarge the
    BH family with tests that cannot be significant, which inflates every q-value in
    that family - conservative, but arbitrarily so.
    """
    drop = set()
    if not (drop_non_bacterial or drop_overview):
        return drop
    with open(path) as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            if drop_non_bacterial and r["is_non_bacterial"] == "yes":
                drop.add(r["pathway_id"])
            if drop_overview and r["is_overview"] == "yes":
                drop.add(r["pathway_id"])
    return drop


def sample_key(col):
    """Sort key giving MF6-, 27D6-, 34F7-, MF6+, 27D6+, 34F7+ ; reps ascending."""
    strain, cult, _day, rep = col.split("_")
    return (0 if cult == "alone" else 1, STRAIN_ORDER.index(strain), int(rep))


def display_name(col):
    """STRAIN_Chlamy_REP / STRAIN_REP - the naming the figure code parses for
    the strain and Chlamy color strips ('_Chlamy_' in the name, strain before
    the first underscore)."""
    strain, cult, _day, rep = col.split("_")
    return f"{strain}_Chlamy_{rep}" if cult == "withC" else f"{strain}_{rep}"





ABSENT_IN_MUTANTS = set()


def contig_of_protein(gff):
    out = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if "locus_tag" in a:
                out[a["locus_tag"]] = f[0]
    return out


def load(args):
    m = pd.read_csv(args.matrix, sep="\t")
    cols = sorted([c for c in m.columns
                   if f"_d{args.day}_" in c and c not in args.exclude],
                  key=sample_key)
    if len(cols) != args.expect:
        sys.exit(f"expected {args.expect} day-{args.day} samples, found "
                 f"{len(cols)}: {cols}")

    X = m.set_index("Protein")[cols].astype(float).replace(0.0, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.log2(X)
    absent = np.asarray(
        L.index.map(lambda p: contig_of_protein(args.gff).get(p) in ABSENT_IN_MUTANTS),
        dtype=bool)
    core = L.notna().all(axis=1) & ~absent
    print(f"excluded           : {int(absent.sum()):,} proteins on "
          f"{'/'.join(sorted(ABSENT_IN_MUTANTS))} (absent from both mutant genomes)")
    print(f"day-{args.day} samples   : {len(cols)}  "
          f"(excluded: {', '.join(args.exclude) or 'none'})")
    print(f"CORE proteins      : {int(core.sum()):,} of {len(L):,} "
          f"(real value in all {len(cols)})")

    offset = L.loc[core].median(axis=0)
    N = L.loc[core].sub(offset, axis=1)
    print(f"median offsets     : {offset.min():.2f} .. {offset.max():.2f} "
          f"(spread {offset.max() - offset.min():.2f} log2)")

    ann = m.set_index("Protein").loc[N.index, ["Gene", "Description"]]
    return N, ann, cols


def cluster(N, k):
    """Per-protein z-score, Ward/euclidean row linkage, maxclust cut into k,
    relabelled 1..k in the order the clusters first appear down the plot."""
    z = N.sub(N.mean(axis=1), axis=0)
    sd = N.std(axis=1).replace(0, np.nan)
    z = z.div(sd, axis=0).fillna(0.0)

    link = linkage(z.values, metric="euclidean", method="ward")
    raw = fcluster(link, t=k, criterion="maxclust")
    leaf_order = leaves_list(link)

    remap = {}
    for leaf in leaf_order:
        c = raw[leaf]
        if c not in remap:
            remap[c] = len(remap) + 1
    clusters = pd.Series([remap[c] for c in raw], index=z.index, name="cluster")

    leaf_rank = pd.Series(np.empty(len(clusters), dtype=int), index=z.index)
    leaf_rank.iloc[leaf_order] = np.arange(len(clusters))
    row_order = (pd.DataFrame({"cluster": clusters, "rank": leaf_rank})
                 .sort_values(["cluster", "rank"]).index)

    sizes = clusters.value_counts().sort_index()
    print(f"row clusters       : k={k} (Ward/euclidean) "
          + ", ".join(f"C{c}:{n}" for c, n in sizes.items()))
    return z, clusters, row_order


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR for a 1-D array of p-values."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def enrich(df, descriptions, exclude=None, min_term_size=1, min_hits=1):
    """df: protein_id, cluster, category (long format -- a protein may appear on
    several rows if it carries several terms). N/n count DISTINCT proteins, not
    protein-term rows, so a multi-term protein never inflates a denominator.
    exclude: categories skipped as tested terms but still counted in the
    background, so N/K/n stay correct.

    min_term_size / min_hits gate what is TESTED, and that matters for more than
    tidiness: BH multiplies each p-value by (family size)/rank, so a term that
    could never reach significance still inflates the q-value of every term it
    shares a cluster with.  On this dataset the untested-but-tested terms were
    519 of 819 KEGG rows and cost 11 real findings - ABC transporters, bacterial
    chemotaxis, the TCA cycle and the two-component system among them.  Terms
    below the floors stay in the background counts (K, N), they simply are not
    tested."""
    exclude = exclude or set()
    N = df["protein_id"].nunique()
    K_t = df.groupby("category")["protein_id"].nunique()
    rows = []
    for c, g in df.groupby("cluster"):
        n = g["protein_id"].nunique()
        counts = g.groupby("category")["protein_id"].nunique()
        recs = []
        for t, k in counts.items():
            if t in exclude:
                continue
            K = int(K_t[t])
            if K < min_term_size or k < min_hits:
                continue
            table = [[k, K - k], [n - k, N - n - (K - k)]]
            _, p = fisher_exact(table, alternative="greater")
            desc, group = descriptions.get(t, (t, ""))
            ids = ";".join(sorted(g.loc[g["category"] == t, "protein_id"]))
            recs.append(dict(cluster=c, category=t, description=desc, group=group,
                             k=k, n=n, K=K, N=N, gene_ratio=k / n,
                             fold_enrichment=(k / n) / (K / N), pvalue=p,
                             protein_ids=ids))
        if not recs:
            continue
        rec = pd.DataFrame(recs)
        rec["FDR"] = bh_fdr(rec["pvalue"].values)
        rows.append(rec)
    cols = ["cluster", "category", "description", "group", "k", "n", "K", "N",
            "gene_ratio", "fold_enrichment", "pvalue", "FDR", "protein_ids"]
    if not rows:
        return pd.DataFrame(columns=cols)
    out = pd.concat(rows, ignore_index=True)
    return out[cols].sort_values(["cluster", "FDR", "pvalue"]).reset_index(drop=True)


def cog_terms(ann_egg):
    """protein_id, category -- one row per protein per COG letter (eggNOG packs
    them into an unseparated string, e.g. 'IQ')."""
    rows = []
    for pid, cats in ann_egg["COG_category"].fillna("").items():
        for letter in str(cats):
            if letter in COG_CATEGORIES:
                rows.append((pid, letter))
    return pd.DataFrame(rows, columns=["protein_id", "category"]).drop_duplicates()


def pathway_terms(ann_egg, drop, known=None):
    """protein_id, category -- one row per protein per KEGG pathway. eggNOG's
    KEGG_Pathway lists each pathway twice, as ko##### and map#####; keep the
    ko-flavoured (organism-agnostic reference) ids only, so nothing is
    double-counted."""
    rows = []
    for pid, paths in ann_egg["KEGG_Pathway"].fillna("").items():
        s = str(paths)
        if not s or s == "-":
            continue
        for t in {t for t in s.split(",") if t.startswith("ko")}:
            if t in drop or (known is not None and t not in known):
                continue
            rows.append((pid, t))
    return pd.DataFrame(rows, columns=["protein_id", "category"]).drop_duplicates()


def run_one(name, terms, descriptions, cl, out_tsv, source, fdr, exclude=None,
            min_term_size=1, min_hits=1):
    df = cl.merge(terms, on="protein_id", how="inner")
    if df.empty:
        sys.exit(f"no {name} terms found for the clustered proteins")
    N = df["protein_id"].nunique()
    n_terms = df.loc[~df["category"].isin(exclude or set()), "category"].nunique()
    print(f"\n[{name}] clustered proteins: {len(cl)} | universe N={N} "
          f"({n_terms} tested terms)")

    res = enrich(df, descriptions, exclude=exclude,
                 min_term_size=min_term_size, min_hits=min_hits)
    header = (
        f"# per-cluster {name} enrichment, source: {source}\n"
        "# test: one-sided Fisher's exact (over-representation), BH-FDR within each cluster\n"
        f"# universe: {N} clustered & {name}-annotated proteins; gene_ratio = k/n; "
        "fold_enrichment = (k/n)/(K/N)\n"
        "# k=cluster proteins with term, n=annotated cluster size, "
        "K=universe proteins with term, N=universe size (all distinct-protein counts;\n"
        "# a protein with >1 term for this term set counts fully toward EVERY term it carries)\n"
        f"# tested only where K >= {min_term_size} and k >= {min_hits}; smaller terms stay in the\n"
        "# background counts but are not tested, so they do not enlarge the BH family\n"
    )
    with open(out_tsv, "w") as fh:
        fh.write(header)
        res.to_csv(fh, sep="\t", index=False)
    print(f"  wrote {out_tsv} ({len(res)} cluster x term tests)")

    sig = res[res["FDR"] <= fdr]
    for c in sorted(res["cluster"].unique()):
        d = sig[sig["cluster"] == c].sort_values("FDR")
        top = (f"top: {d.iloc[0]['category']} ({d.iloc[0]['description'][:32]}, "
               f"FDR={d.iloc[0]['FDR']:.1e})" if not d.empty else "no significant terms")
        print(f"  C{c}: {len(d)} sig terms | {top}")
    return res


def main() -> int:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        matrix=Path(sm.input.matrix),
        eggnog=Path(sm.input.eggnog),
        gff=Path(sm.input.gff),
        kegg_names=Path(sm.input.kegg),
        kegg_classes=Path(sm.input.kegg_classes),
        outdir=Path(sm.params.outdir),
        absent_replicon=sm.params.absent,
        prefix="core",
        day=int(sm.params.day),
        exclude=list(sm.params.exclude),
        expect=int(sm.params.expect),
        clusters=int(sm.params.k),
        fdr=float(sm.params.fdr),
        drop_retired=bool(sm.params.drop_ret),
        drop_non_bacterial=bool(sm.params.drop_nb),
        drop_kegg_overview=bool(sm.params.drop_ov),
        min_term_size=int(sm.params.min_term),
        min_hits=int(sm.params.min_hits),
    )
    global ABSENT_IN_MUTANTS
    ABSENT_IN_MUTANTS = {args.absent_replicon}
    args.outdir.mkdir(parents=True, exist_ok=True)

    N, ann, cols = load(args)
    z, clusters, row_order = cluster(N, args.clusters)

    egg = pd.read_csv(args.eggnog, sep="\t", skiprows=4)
    egg.columns = [c.lstrip("#") for c in egg.columns]
    egg["pid"] = egg["query"].astype(str).str.split("|").str[-1]
    egg = egg.drop_duplicates("pid").set_index("pid")
    keep = ["COG_category", "KEGG_ko", "KEGG_Pathway", "Preferred_name"]
    ann_egg = egg.reindex(N.index)[keep]
    print(f"eggNOG annotated   : {int(ann_egg['COG_category'].notna().sum()):,} "
          f"of {len(N):,} core proteins have a COG category")

    table = pd.concat([ann, ann_egg], axis=1).loc[row_order].copy()
    table.insert(0, "protein_id", table.index)
    table["cluster"] = clusters.loc[row_order].values
    for orig in cols:
        table[display_name(orig)] = N.loc[row_order, orig].values

    clusters_tsv = args.outdir / f"{args.prefix}_clusters.tsv"
    header = (
        f"# shared/complete-case day-{args.day} proteome, {len(N)} proteins "
        f"(real value in all {len(cols)} samples)\n"
        f"# excluded runs: {', '.join(args.exclude) or 'none'}\n"
        f"# clustering: rows ward/euclidean on the per-protein z-score, cut to "
        f"k={args.clusters}; columns fixed order (not clustered): "
        "MF6-,27D6-,34F7-,MF6+,27D6+,34F7+, rep ascending\n"
        "# sample columns renamed STRAIN_Chlamy_REP (+Chlamy) / STRAIN_REP (-Chlamy); "
        "values are the median-centred log2 MaxLFQ abundances\n"
    )
    with open(clusters_tsv, "w") as fh:
        fh.write(header)
        table.to_csv(fh, sep="\t", index=False)
    print(f"\nwrote {clusters_tsv} ({len(table)} rows)")
    print("column order:", " ".join(display_name(c) for c in cols))

    cl = table[["protein_id", "cluster"]]
    run_one("COG", cog_terms(ann_egg), COG_CATEGORIES, cl,
            args.outdir / f"{args.prefix}_cog_enrichment.tsv", clusters_tsv,
            args.fdr, exclude=None,
            min_term_size=args.min_term_size, min_hits=args.min_hits)

    kegg_desc = {}
    if args.kegg_names.is_file():
        nm = pd.read_csv(args.kegg_names, sep="\t")
        kegg_desc = {pid: (name, "") for pid, name in zip(nm["pathway_id"], nm["name"])}
    else:
        print(f"  WARNING: {args.kegg_names} missing; pathway labels stay bare ko ids")
    drop = (kegg_drop_set(args.kegg_classes, args.drop_non_bacterial,
                          args.drop_kegg_overview) if args.kegg_classes else set())
    print(f"KEGG filter        : {len(drop)} maps dropped by class "
          f"({'non-bacterial ' if args.drop_non_bacterial else ''}"
          f"{'overview' if args.drop_kegg_overview else ''})")
    known = (kegg_known_set(args.kegg_classes)
             if args.kegg_classes and args.drop_retired else None)
    if known is not None:
        print(f"                     + retired maps dropped (not among "
              f"{len(known)} current KEGG maps)")
    run_one("KEGG-pathway", pathway_terms(ann_egg, drop, known),
            kegg_desc, cl, args.outdir / f"{args.prefix}_pathway_enrichment.tsv",
            clusters_tsv, args.fdr,
            min_term_size=args.min_term_size, min_hits=args.min_hits)
    return 0


if __name__ == "__main__":
    sys.exit(main())
