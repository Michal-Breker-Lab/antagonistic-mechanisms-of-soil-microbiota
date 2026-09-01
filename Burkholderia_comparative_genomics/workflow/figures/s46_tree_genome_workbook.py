#!/usr/bin/env python3
"""Build the Burkholderia tree genome + metadata + tree-statistics workbook."""
import json, re, os, sys
from pathlib import Path
import pandas as pd
from Bio import Phylo
from io import StringIO

P = Path(os.path.expanduser("~/Moshe/Projects/Burkholderia_c3_Pangenome"))
R = P / "rebuild"
T = R / "tables"
RES = R / "results"
OUT = P / "tables" / "Burkholderia_tree_genomes_and_metadata.xlsx"

def tsv(p, **kw):
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, **kw)

# ---------- trees ----------
def load_tree(path):
    return Phylo.read(str(path), "newick")

chr1 = load_tree(T / "tree_chr1.treefile")
pc3  = load_tree(RES / "tree_pc3.treefile")

def tip_table(tree, label):
    rows = []
    # depth from root in branch-length units
    depths = tree.depths(unit_branch_lengths=False)
    ndepth = tree.depths(unit_branch_lengths=True)
    for t in tree.get_terminals():
        rows.append({"accession": t.name,
                     f"{label}_tip_branch_length": t.branch_length,
                     f"{label}_root_to_tip_distance": depths.get(t),
                     f"{label}_nodes_from_root": ndepth.get(t)})
    return pd.DataFrame(rows)

chr1_tips = tip_table(chr1, "chr1")
pc3_tips  = tip_table(pc3, "pC3")

# ---------- display rooting: reuse the project's own conventions ----------
sys.path.insert(0, str(R / "scripts"))
import tree_display as td

_c3raw = tsv(RES / "c3_calls_all_genomes.tsv")
c3_present_map = dict(zip(_c3raw["accession"],
                          _c3raw["c3_present"].str.lower().eq("true")))
_cen = tsv(T / "replicon_census.tsv")
org_map = dict(zip(_cen["accession"], _cen["organism_name"]))

chr1_disp = load_tree(T / "tree_chr1.treefile")
chr1_disp.root_at_midpoint()
chr1_disp.ladderize()
ingroup, outgroup = td.ingroup_root(chr1_disp, c3_present_map, org_map,
                                    verbose=False)
# `ingroup` is stale after root_with_outgroup re-roots in place (fig1 uses only
# `outgroup` for the same reason), so define the ingroup as its complement.
out_names = {t.name for t in outgroup.get_terminals()}
in_names = {t.name for t in chr1_disp.get_terminals()} - out_names
nearclones = td.nearclone_clades(chr1_disp, c3_present_map, org_map,
                                 skip=[outgroup])
wedge = {}
for i, cl in enumerate(nearclones, 1):
    lbl = f"wedge{i:02d} ({len(cl.get_terminals())} tips)"
    for t in cl.get_terminals():
        wedge[t.name] = lbl

_depth = chr1_disp.depths(unit_branch_lengths=False)
root_dist = {t.name: _depth.get(t) for t in chr1_disp.get_terminals()}
chr1_tips["chr1_root_to_tip_distance"] = chr1_tips["accession"].map(root_dist)
chr1_tips["chr1_tree_position"] = chr1_tips["accession"].map(
    lambda a: "ingroup" if a in in_names else "outgroup")
chr1_tips["chr1_display_collapsed_wedge"] = chr1_tips["accession"].map(
    lambda a: wedge.get(a, ""))

# ---------- per-genome metadata ----------
census = tsv(T / "replicon_census.tsv")
hosts  = tsv(T / "host_categories.tsv")
qc     = tsv(RES / "annotation_qc.tsv")
c3     = tsv(RES / "c3_calls_all_genomes.tsv")
clone  = tsv(T / "clone_cluster.tsv")
gs     = tsv(R / "genus_species.tsv", header=None,
             names=["accession", "genus_parsed", "species_parsed"])
gapf   = tsv(T / "chr1_alignment_gap_fraction.tsv", header=None,
             names=["accession", "chr1_core_gap_fraction"])

# ANI to MF6 (best value per genome)
ani = tsv(T / "MF6_ani_raw.tsv")
ani["accession"] = ani["Ref_file"].str.replace(r"^genomes/", "", regex=True) \
                                 .str.replace(r"\.fna$", "", regex=True)
ani["ANI"] = pd.to_numeric(ani["ANI"], errors="coerce")
ani_best = (ani.sort_values("ANI", ascending=False)
              .drop_duplicates("accession")[["accession", "ANI",
                                             "Align_fraction_ref",
                                             "Align_fraction_query"]]
              .rename(columns={"ANI": "ani_to_MF6",
                               "Align_fraction_ref": "MF6_align_frac_ref",
                               "Align_fraction_query": "MF6_align_frac_query"}))

# assembly level where the bac120 harvest recorded it
meta = json.load(open(P / "bac120" / "genome_metadata.json"))
lvl = pd.DataFrame([{"accession": k, "assembly_level": v.get("level", "")}
                    for k, v in meta.items()])

# toxin / RHS carriage per genome (tier-1 loci carried, out of 12)
tox = tsv(T / "toxin12_carriers.tsv")
tox["tier1_carrier"] = tox["tier1_carrier"].str.lower().eq("true")
tox_n = (tox.groupby("accession")["tier1_carrier"].sum()
            .rename("n_toxin12_loci_tier1").reset_index())
rhs = tsv(T / "rhs_search_per_genome.tsv")

# ---------- assemble ----------
genomes = [x.strip() for x in open(R / "genome_list_full.txt") if x.strip()]
df = pd.DataFrame({"accession": genomes})

df = (df.merge(census, on="accession", how="left")
        .merge(gs, on="accession", how="left")
        .merge(hosts.drop(columns=["organism_name"]), on="accession", how="left")
        .merge(lvl, on="accession", how="left")
        .merge(qc, on="accession", how="left")
        .merge(c3, on="accession", how="left")
        .merge(clone, on="accession", how="left")
        .merge(ani_best, on="accession", how="left")
        .merge(gapf, on="accession", how="left")
        .merge(tox_n, on="accession", how="left")
        .merge(chr1_tips, on="accession", how="left")
        .merge(pc3_tips, on="accession", how="left"))

df["on_chr1_tree"] = df["accession"].isin(set(chr1_tips["accession"]))
df["on_pC3_tree"]  = df["accession"].isin(set(pc3_tips["accession"]))

# exclusion reasons (D13 / D20 / QC)
d13 = {"GCA_040954445.1": "99.18", "GCF_050430075.1": "95.97",
       "GCF_034047095.1": "95.97", "GCF_039852015.1": "95.62",
       "GCF_004842085.1": "93.97", "MF7": "93.96",
       "GCF_003812585.1": "93.86", "GCF_003568605.1": "91.06"}
d20 = {"GCA_059696275.1", "GCF_050955445.1"}

def reason(a, on):
    if on:
        return ""
    if a in d13:
        return f"D13: >90% gaps in chr1 core alignment ({d13[a]}%)"
    if a in d20:
        return "D20: holds 0 of 1,018 chr1 alignment families; never entered the alignment"
    return "not in chr1 core alignment"

df["chr1_tree_exclusion_reason"] = [reason(a, o) for a, o
                                    in zip(df["accession"], df["on_chr1_tree"])]

cols = ["accession", "organism_name", "genus_parsed", "species_parsed", "strain",
        "assembly_level", "host_category", "raw_host", "raw_isolation_source",
        "raw_host_disease", "geo_loc_name", "evidence_x" if "evidence_x" in df else "evidence",
        "total_length", "gc_percent", "n_replicons_total", "largest", "second",
        "third", "architecture", "ncbi_labelled_chromosomes",
        "ncbi_labelled_plasmids", "n_large_200kb", "n_large_500kb",
        "sequencing_tech", "checkm_completeness", "checkm_contamination",
        "qc_pass", "qc_reason", "CDS", "tRNA", "rRNA", "contigs",
        "c3_present", "n_secondary_large", "cluster_id", "cluster_size",
        "is_representative", "ani_to_MF6", "MF6_align_frac_ref",
        "n_toxin12_loci_tier1", "chr1_core_gap_fraction",
        "on_chr1_tree", "chr1_tree_position", "chr1_tip_branch_length",
        "chr1_root_to_tip_distance", "chr1_display_collapsed_wedge",
        "chr1_tree_exclusion_reason", "on_pC3_tree", "pC3_tip_branch_length",
        "pC3_root_to_tip_distance"]
cols = [c for c in cols if c in df.columns]
df = df[cols]

rename = {"genus_parsed": "genus", "species_parsed": "species",
          "evidence_x": "host_evidence", "evidence": "host_evidence",
          "cluster_id": "clone_cluster_id", "cluster_size": "clone_cluster_size",
          "is_representative": "clone_cluster_representative",
          "largest": "replicon1_bp", "second": "replicon2_bp",
          "third": "replicon3_bp"}
df = df.rename(columns=rename)

# numeric coercion
num = ["total_length", "gc_percent", "n_replicons_total", "replicon1_bp",
       "replicon2_bp", "replicon3_bp", "ncbi_labelled_chromosomes",
       "ncbi_labelled_plasmids", "n_large_200kb", "n_large_500kb",
       "checkm_completeness", "checkm_contamination", "CDS", "tRNA", "rRNA",
       "contigs", "n_secondary_large", "clone_cluster_id", "clone_cluster_size",
       "ani_to_MF6", "MF6_align_frac_ref", "n_toxin12_loci_tier1",
       "chr1_core_gap_fraction", "chr1_tip_branch_length",
       "chr1_root_to_tip_distance", "pC3_tip_branch_length",
       "pC3_root_to_tip_distance"]
for c in num:
    if c in df:
        df[c] = pd.to_numeric(df[c], errors="coerce")
for c in ["qc_pass", "c3_present", "clone_cluster_representative"]:
    if c in df:
        df[c] = df[c].astype(str).str.title()

df = df.sort_values("accession").reset_index(drop=True)
on_tree = df[df["on_chr1_tree"]].copy()
off_tree = df[~df["on_chr1_tree"]].copy()
on_tree = on_tree.drop(columns=["on_chr1_tree", "chr1_tree_exclusion_reason"])
off_tree = off_tree.drop(columns=[c for c in off_tree.columns
                                  if c.startswith(("chr1_tip", "chr1_root",
                                                   "chr1_tree_position",
                                                   "chr1_display"))
                                  or c == "on_chr1_tree"])

# ---------- node support ----------
def support_table(tree, label):
    rows = []
    for i, cl in enumerate(tree.get_nonterminals(), 1):
        s = cl.confidence
        if s is None and cl.name:
            try:
                s = float(cl.name)
            except ValueError:
                s = None
        rows.append({"node": f"{label}_N{i}", "ufboot_support": s,
                     "branch_length": cl.branch_length,
                     "n_descendant_tips": len(cl.get_terminals()),
                     "is_root": cl is tree.root})
    return pd.DataFrame(rows)

sup1 = support_table(chr1, "chr1")
sup2 = support_table(pc3, "pC3")

def support_summary(s, label):
    v = s["ufboot_support"].dropna()
    bl = s["branch_length"].dropna()
    return pd.DataFrame({
        "tree": label,
        "metric": ["internal nodes", "nodes with support value",
                   "median UFBoot", "mean UFBoot",
                   "nodes UFBoot >= 95 (n)", "nodes UFBoot >= 95 (%)",
                   "nodes UFBoot >= 80 (n)", "nodes UFBoot >= 80 (%)",
                   "nodes UFBoot < 50 (n)",
                   "internal branches with length 0 (n)",
                   "median internal branch length"],
        "value": [len(s), len(v), v.median(), round(v.mean(), 2),
                  int((v >= 95).sum()), round(100 * (v >= 95).mean(), 1),
                  int((v >= 80).sum()), round(100 * (v >= 80).mean(), 1),
                  int((v < 50).sum()), int((bl == 0).sum()),
                  round(bl.median(), 8)]})

# ---------- bac120 marker tree ----------
B = P / "bac120"
bac_tree = load_tree(B / "tables" / "iq_bac120.treefile")
bac_tips = tip_table(bac_tree, "bac120")

bmeta = json.load(open(B / "genome_metadata.json"))
bdf = pd.DataFrame([{"accession": k,
                     "organism_name": v.get("organism", ""),
                     "assembly_level": v.get("level", ""),
                     "raw_host": v.get("host", ""),
                     "raw_isolation_source": v.get("isolation_source", ""),
                     "env_medium": v.get("env_medium", "")}
                    for k, v in bmeta.items()])

tn = tsv(B / "tables" / "lab_derivatives_excluded.tsv")
tn = tn.rename(columns={"strain": "tn_strain", "reason": "lab_derivative_reason",
                        "assembly_level": "_lvl", "n_contigs": "_nc"})
tn["is_lab_derivative"] = True
tn = tn[["accession", "is_lab_derivative", "parent_accession", "parent_strain",
         "bioproject", "bioproject_title", "lab_derivative_reason"]]

bani = tsv(B / "tables" / "mf6_ani_full.tsv")
bani["accession"] = bani["Ref_file"].str.replace(r"^genomes/", "", regex=True) \
                                    .str.replace(r"\.fna$", "", regex=True)
bani["ANI"] = pd.to_numeric(bani["ANI"], errors="coerce")
bani_best = (bani.sort_values("ANI", ascending=False).drop_duplicates("accession")
               [["accession", "ANI"]].rename(columns={"ANI": "ani_to_MF6"}))

bc3 = tsv(B / "tables" / "c3_calls_new.tsv")[["accession", "c3_present", "evidence"]] \
        .rename(columns={"evidence": "c3_evidence"})
old_c3 = tsv(RES / "c3_calls_all_genomes.tsv")[["accession", "c3_present", "evidence"]] \
           .rename(columns={"evidence": "c3_evidence"})
allc3 = pd.concat([old_c3, bc3]).drop_duplicates("accession", keep="last")

bac_df = (bac_tips.merge(bdf, on="accession", how="left")
                  .merge(tn, on="accession", how="left")
                  .merge(bani_best, on="accession", how="left")
                  .merge(allc3, on="accession", how="left"))
bac_df["is_lab_derivative"] = bac_df["is_lab_derivative"].fillna(False)
bac_df["also_on_chr1_tree"] = bac_df["accession"].isin(set(chr1_tips["accession"]))
bac_df["in_773_genome_set"] = bac_df["accession"].isin(set(genomes))
bac_df["ani_to_MF6"] = pd.to_numeric(bac_df["ani_to_MF6"], errors="coerce")
bac_df = bac_df[["accession", "organism_name", "assembly_level", "raw_host",
                 "raw_isolation_source", "env_medium", "c3_present", "c3_evidence",
                 "ani_to_MF6", "is_lab_derivative", "parent_accession",
                 "parent_strain", "bioproject", "bioproject_title",
                 "lab_derivative_reason", "also_on_chr1_tree", "in_773_genome_set",
                 "bac120_tip_branch_length", "bac120_root_to_tip_distance",
                 "bac120_nodes_from_root"]].sort_values("accession").reset_index(drop=True)

sup3 = support_table(bac_tree, "bac120")

# ---------- tree statistics ----------
def iqtree_stats(path):
    txt = Path(path).read_text()
    g = lambda p: (re.search(p, txt).group(1).strip() if re.search(p, txt) else "")
    return {
        "IQ-TREE version": g(r"^(IQ-TREE [\d.]+)"),
        "alignment": g(r"Input file name: (.+)"),
        "taxa": g(r"Input data: (\d+) sequences"),
        "alignment sites": g(r"sequences with (\d+ (?:nucleotide|amino-acid) sites)"),
        "constant sites": g(r"Number of constant sites: (\d+ \(= [\d.]+% of all sites\))"),
        "parsimony-informative sites": g(r"Number of parsimony informative sites: (\d+)"),
        "distinct site patterns": g(r"Number of distinct site patterns: (\d+)"),
        "substitution model": g(r"Model of substitution: (\S+)"),
        "rate categories": g(r"Site proportion and rates:\s+(.+)"),
        "log-likelihood": g(r"Log-likelihood of the tree: (-[\d.]+ \(s\.e\. [\d.]+\))"),
        "unconstrained logL": g(r"Unconstrained log-likelihood \(without tree\): (-[\d.]+)"),
        "AIC": g(r"Akaike information criterion \(AIC\) score: ([\d.]+)"),
        "BIC": g(r"Bayesian information criterion \(BIC\) score: ([\d.]+)"),
        "total tree length": g(r"Total tree length \(sum of branch lengths\): ([\d.]+)"),
        "sum of internal branches": g(r"Sum of internal branch lengths: (.+)"),
        "near-zero internal branches": g(r"WARNING: (\d+) near-zero internal branches"),
        "branch support": "ultrafast bootstrap (UFBoot), 1000 replicates",
    }

st1 = iqtree_stats(T / "tree_chr1.iqtree")
st2 = iqtree_stats(RES / "tree_pc3.iqtree")
st3 = iqtree_stats(P / "bac120" / "tables" / "iq_bac120.iqtree")
stats = pd.DataFrame({"statistic": list(st1.keys()),
                      "chr1 core species tree": list(st1.values()),
                      "pC3 core tree": [st2.get(k, "") for k in st1],
                      "bac120 marker tree": [st3.get(k, "") for k in st1]})

extra = pd.DataFrame({
    "statistic": ["tips on tree", "genomes in full set", "genomes not on this tree",
                  "UFBoot convergence (bootstrap correlation)",
                  "pangenome software", "core threshold",
                  "core genes (99-100% of strains)", "soft-core genes (95-99%)",
                  "shell genes (15-95%)", "cloud genes (<15%)", "total gene clusters",
                  "tree file", "IQ-TREE report"],
    "chr1 core species tree": [
        len(chr1_tips), 773, 10, "0.983 (accepted; --bcor 0.98)",
        "Panaroo (strict)", "strict",
        "65", "953", "2947", "142413", "146378",
        "rebuild/tables/tree_chr1.treefile",
        "rebuild/tables/tree_chr1.iqtree"],
    "pC3 core tree": [
        len(pc3_tips), 773, "520 (pC3-negative; the tree covers pC3 carriers only)", "",
        "Panaroo (moderate)", "moderate",
        "24", "31", "1347", "18628", "20030",
        "rebuild/results/tree_pc3.treefile",
        "rebuild/results/tree_pc3.iqtree"],
    "bac120 marker tree": [
        len(bac_tips), 791, "n/a (separate NCBI harvest, not the 773-genome set)", "",
        "GTDB-Tk bac120 marker set", "120 concatenated bacterial marker proteins",
        "n/a", "n/a", "n/a", "n/a", "n/a",
        "bac120/tables/iq_bac120.treefile",
        "bac120/tables/iq_bac120.iqtree"]})
stats = pd.concat([stats, extra], ignore_index=True)

# ---------- composition ----------
sp = (on_tree.groupby("organism_name").size().rename("n_tips")
        .reset_index().sort_values("n_tips", ascending=False))
sp["percent_of_tips"] = (100 * sp["n_tips"] / len(on_tree)).round(2)

hc = (on_tree.groupby("host_category").size().rename("n_tips").reset_index()
        .sort_values("n_tips", ascending=False))
hc["percent_of_tips"] = (100 * hc["n_tips"] / len(on_tree)).round(2)
hc_c3 = on_tree.pivot_table(index="host_category", columns="c3_present",
                            values="accession", aggfunc="count").fillna(0).astype(int)
hc = hc.merge(hc_c3.reset_index(), on="host_category", how="left")

gen = (on_tree.groupby("genus").size().rename("n_tips").reset_index()
         .sort_values("n_tips", ascending=False))
gen["percent_of_tips"] = (100 * gen["n_tips"] / len(on_tree)).round(2)

arch = (on_tree.groupby("architecture").size().rename("n_tips")
          .reset_index().sort_values("n_tips", ascending=False))

# ---------- readme ----------
readme = pd.DataFrame({
"Sheet": ["README", "Genomes_on_chr1_tree", "Genomes_excluded",
          "Tree_statistics", "Node_support_summary", "chr1_node_support",
          "pC3_node_support", "Genomes_bac120_tree", "bac120_node_support",
          "Species_composition", "Genus_composition",
          "Host_composition", "Replicon_architecture"],
"Contents": [
 "This sheet: provenance and column notes.",
 f"{len(on_tree)} genomes that are tips of the chromosome-1 core species tree, with full metadata, "
 "assembly QC, host/isolation metadata, replicon architecture, pC3 (c3) call, clone cluster, "
 "ANI to MF6, chr1 core-alignment gap fraction, and tip branch lengths.",
 f"{len(off_tree)} genomes in the 773-genome set that are NOT tips of the chr1 tree, each with the reason.",
 "Per-tree phylogenetic statistics: alignment, model, likelihood, tree length, support scheme, pangenome core.",
 "Branch-support distributions for both trees.",
 "Every internal node of the chr1 tree: UFBoot support, branch length, descendant tip count.",
 "Every internal node of the pC3 tree: UFBoot support, branch length, descendant tip count.",
 "791 genomes on the complementary GTDB-Tk bac120 marker tree (a separate, later NCBI "
 "harvest — NOT the same genome set as the chr1 tree). Flags the 37 transposon mutants "
 "of one B. sola strain (PRJEB40633) that inherit the parent's isolate metadata.",
 "Every internal node of the bac120 tree: UFBoot support, branch length, descendant tip count.",
 "Tip counts per organism_name on the chr1 tree (sampling-bias check).",
 "Tip counts per genus on the chr1 tree.",
 "Tip counts per host category on the chr1 tree, split by pC3 presence.",
 "Tip counts per replicon architecture class on the chr1 tree."]})

prov = pd.DataFrame({"Item": [
 "Project", "Source directory", "Tree (primary)", "Tree (pC3)", "Tree (bac120)",
 "Genome set", "Dereplication", "Metadata sources", "Generated"],
 "Value": [
 "Burkholderia_c3_Pangenome — non-dereplicated rebuild",
 str(P),
 "rebuild/tables/tree_chr1.treefile  (= tables/chr1_core.treefile, identical) — 763 tips",
 "rebuild/results/tree_pc3.treefile — 253 tips",
 "bac120/tables/iq_bac120.treefile — 791 tips (complementary marker backbone, different genome set)",
 "773 = all complete Burkholderia sensu lato genomes + lab isolates MF6 and MF7",
 "None. Clone clusters (99% ANI) are carried as a covariate (clone_cluster_id), not a filter.",
 "replicon_census.tsv, host_categories.tsv, annotation_qc.tsv, c3_calls_all_genomes.tsv, "
 "clone_cluster.tsv, genus_species.tsv, chr1_alignment_gap_fraction.tsv, MF6_ani_raw.tsv, "
 "toxin12_carriers.tsv, bac120/genome_metadata.json (assembly level), tree_chr1.iqtree, tree_pc3.iqtree",
 pd.Timestamp.today().strftime("%Y-%m-%d")]})

notes = pd.DataFrame({"Column": [
 "c3_present", "clone_cluster_id", "ani_to_MF6", "chr1_core_gap_fraction",
 "assembly_level", "checkm_completeness", "n_toxin12_loci_tier1",
 "chr1_tip_branch_length", "chr1_root_to_tip_distance", "chr1_tree_position",
 "chr1_display_collapsed_wedge", "host_category",
 "bac120 UFBoot support", "is_lab_derivative"],
 "Note": [
 "pC3 (third-replicon) presence call; 253 of 773 positive.",
 "99% ANI clone cluster; cluster_size 1 means the genome is unique in the set.",
 "Best ANI of this genome against MF6 (skani). Blank = below the search reporting threshold.",
 "Fraction of the chromosome-1 core alignment that is gap in this genome. "
 "D13 excluded every genome above 0.90.",
 "From the bac120 NCBI harvest; blank where that harvest did not cover the accession. "
 "The whole set is complete/chromosome-level by construction.",
 "CheckM completeness where NCBI reported it; blank is common for complete genomes.",
 "Number of the 12 MF6 toxin/immunity loci carried at tier-1 stringency.",
 "Terminal branch length in substitutions/site from IQ-TREE.",
 "IQ-TREE writes an UNROOTED tree. This distance is measured on the tree after the "
 "project's display rooting (midpoint, then re-rooted on the ingroup stem, as in "
 "tree_display.ingroup_root used by Figures 1/4/10/11). It is a display quantity, "
 "not a molecular clock.",
 "Ingroup = the largest clade containing every pC3 carrier and no non-Burkholderia-"
 "labelled tip (638 tips, 100% UFBoot). Defined topologically because the organism "
 "labels are polyphyletic here.",
 "Near-clone wedge this tip is collapsed into for figure display (blank = drawn "
 "individually). Collapsing is display-only and never hides a pC3 carrier.",
 "Curated from raw_host / raw_isolation_source / raw_host_disease; see host_evidence.",
 "bac120 support is weak by design of the marker set: only 735 of 5,010 sites are "
 "parsimony-informative, so median UFBoot is 37 and 532 internal branches are near-zero. "
 "It resolves deep genus-level structure, not within-species relationships — read the "
 "chr1 core tree for those.",
 "37 of the 791 bac120 tips are isogenic transposon mutants of one B. sola strain "
 "(PRJEB40633) that inherit the parent BioSample's metadata. They inflate any per-genome "
 "count; drop them for display, keep the parent GCF_905400185.1."]})

# ---------- write ----------
OUT.parent.mkdir(parents=True, exist_ok=True)
with pd.ExcelWriter(OUT, engine="xlsxwriter") as xw:
    readme.to_excel(xw, sheet_name="README", index=False, startrow=0)
    prov.to_excel(xw, sheet_name="README", index=False, startrow=len(readme) + 3)
    notes.to_excel(xw, sheet_name="README", index=False, startrow=len(readme) + len(prov) + 6)
    on_tree.to_excel(xw, sheet_name="Genomes_on_chr1_tree", index=False)
    off_tree.to_excel(xw, sheet_name="Genomes_excluded", index=False)
    stats.to_excel(xw, sheet_name="Tree_statistics", index=False)
    pd.concat([support_summary(sup1, "chr1 core species tree"),
               support_summary(sup2, "pC3 core tree"),
               support_summary(sup3, "bac120 marker tree")], ignore_index=True) \
      .to_excel(xw, sheet_name="Node_support_summary", index=False)
    sup1.to_excel(xw, sheet_name="chr1_node_support", index=False)
    sup2.to_excel(xw, sheet_name="pC3_node_support", index=False)
    bac_df.to_excel(xw, sheet_name="Genomes_bac120_tree", index=False)
    sup3.to_excel(xw, sheet_name="bac120_node_support", index=False)
    sp.to_excel(xw, sheet_name="Species_composition", index=False)
    gen.to_excel(xw, sheet_name="Genus_composition", index=False)
    hc.to_excel(xw, sheet_name="Host_composition", index=False)
    arch.to_excel(xw, sheet_name="Replicon_architecture", index=False)

    wb = xw.book
    hdr = wb.add_format({"bold": True, "bg_color": "#DCE6F1", "border": 1,
                         "text_wrap": True, "valign": "top"})
    wrap = wb.add_format({"text_wrap": True, "valign": "top"})
    for name, frame in [("Genomes_on_chr1_tree", on_tree),
                        ("Genomes_excluded", off_tree),
                        ("Tree_statistics", stats),
                        ("Node_support_summary", None),
                        ("chr1_node_support", sup1),
                        ("pC3_node_support", sup2),
                        ("Genomes_bac120_tree", bac_df),
                        ("bac120_node_support", sup3),
                        ("Species_composition", sp),
                        ("Genus_composition", gen),
                        ("Host_composition", hc),
                        ("Replicon_architecture", arch)]:
        ws = xw.sheets[name]
        if frame is not None:
            for j, c in enumerate(frame.columns):
                ws.write(0, j, c, hdr)
                w = max(len(str(c)), *(frame[c].astype(str).str.len().max(),)) if len(frame) else len(str(c))
                ws.set_column(j, j, min(max(10, w + 2), 46))
            ws.freeze_panes(1, 1)
            ws.autofilter(0, 0, len(frame), len(frame.columns) - 1)
    ws = xw.sheets["README"]
    ws.set_column(0, 0, 30, wrap)
    ws.set_column(1, 1, 110, wrap)

print("wrote", OUT)
print("on-tree:", len(on_tree), "| excluded:", len(off_tree))
print(df["on_chr1_tree"].value_counts().to_dict())
