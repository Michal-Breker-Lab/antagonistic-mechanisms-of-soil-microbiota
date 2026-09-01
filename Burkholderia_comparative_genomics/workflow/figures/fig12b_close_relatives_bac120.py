#!/usr/bin/env python3
"""Figure 12b - MF6 and its closest relatives, pruned out of the bac120 MARKER tree.

This is the rebuild Moshe asked for: "so MF7 can be in the tree". It is a prune of
bac120/tables/<the bac120 tree>, the non-dereplicated 791-genome GTDB-Tk marker
phylogeny, NOT of chr1_core.treefile.

Selection is the clade N splits above MF6 (--splits, default 6), and the 37
isogenic transposon mutants of B. sola R-12632 that BioProject PRJEB40633 put on
that tree are dropped in favour of their parent -- see
tree_display.drop_lab_derivatives.

WHAT THE REBUILD FIXES, and it is worth being precise because it is exactly two
things:
  1. MF7 is a REAL INFERRED TIP. In fig12 it had to be grafted at zero divergence
     with an open star, because its 63-contig draft was dropped from the
     chromosome-1 alignment (D13, 93.96% gaps). On marker genes it is complete:
     118 unique + 2 multi = all 120 bac120 markers, byte-identical counts to MF6.
  2. MF6's two closest relatives are BACK. GCF_016899425.1 (MS389, 99.71% ANI,
     isolated from soybean) and GCF_053209605.1 (B. sola, sputum) were removed by
     dereplication before the chr1 tree was built, so fig12 could only name them
     in a footnote. This set does not dereplicate; both are ordinary tips.

WHAT THE REBUILD COSTS, and this is the honest headline: **the marker tree cannot
resolve MF6 from its closest relatives.** Measured directly on the MSA, the 36
selected genomes carry only 20 DISTINCT bac120 sequences, and the largest
collapse is precisely the group this figure is about:

    MF6, MF7, GCF_016899425.1 (99.71% ANI), GCF_905400185.1 (99.50% ANI)
        -- byte-identical across all 5,010 masked amino-acid columns.

That is not a bug in the run. 120 marker genes at 5,010 aligned positions simply
do not carry enough signal to separate genomes that are >99.5% identical, whereas
the chromosome-1 alignment had 1,337,868 bp to work with. The consequence for
this figure is that the phylogram degenerates to a polytomy exactly where fig12
showed structure. It is drawn honestly rather than resolved arbitrarily, and the
caption says so. If resolution among MF6's neighbourhood is what is wanted, the
instrument is a dedicated core-gene alignment of these ~36 genomes (as in 6.3),
not this tree.

fig12_close_relatives.py is deliberately LEFT IN PLACE and not overwritten: it is
the higher-resolution view of the same clade, and destroying it to gain MF7 would
be a bad trade. The two are complementary.

Inputs: bac120/tables/{iq_bac120.contree | tree_bac120.nwk}, bac120/final_accessions.json,
        rebuild/tables/{MF6_ani_raw.tsv, host_categories.tsv, c3_calls_all_genomes.tsv,
        toxin12_carriers.tsv, toxin12_locus_summary.tsv}
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import tree_display as td  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
BAC = D.parent / "bac120"
FIG.mkdir(exist_ok=True)

# ------------------------------------------------------------------ selection
# SELECTION IS NOW TOPOLOGICAL, not ANI-based. Omri asked (2026-08-25, pointing
# at Figure 13) that this panel show "all genomes from 4 splits above MF6", i.e.
# every descendant of an ancestor node of MF6 rather than everything above an ANI
# cut. The two are not the same question: an ANI cut asks "how similar to MF6",
# a clade asks "what is MF6 nested inside", and only the second is answerable
# from the tree that is being drawn.
#
# Counting matters, and the count that a reader makes off Figure 13 is not the
# topological one, because MF6's first two ancestors sit on branches at IQ-TREE's
# 2.55e-06 floor and are invisible at that scale. Ancestors of MF6, measured:
#     1 up ->  2 tips     5 up ->  49 tips (>=94.28% ANI)
#     2 up ->  5 tips     6 up ->  84 tips (>=93.99% ANI)  <- the node arrowed
#     3 up -> 43 tips     7 up ->  85 tips
#     4 up -> 45 tips     8 up -> 146 tips
# The arrow in the marked-up Figure 13 lands on the 84-tip node (its drawn
# position is r=0.414, angle=31.9 deg; the arrowhead resolves to r=0.45,
# angle=34 deg). The literal 4th topological ancestor is the 45-tip clade, which
# is EXACTLY what the previous version of this figure already showed (>=95% ANI
# selected 45 genomes and they were exactly monophyletic) -- so reading "4
# splits" topologically would be a no-op, and the intended change is the arrowed
# node. Set SPLITS_UP=4 to reproduce the previous panel.
SPLITS_UP = int(sys.argv[sys.argv.index("--splits") + 1]) if "--splits" in sys.argv else 6

# kept only for the ANI column and for reporting; no longer a selection rule
ANI_MIN = 95.0

PROT_COLOR = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
              "#DDCC77", "#CC6677", "#882255", "#AA4499", "#EE7733"]
RING_LOCI = ["MF6_001079", "MF6_002734", "MF6_003184",
             "MF6_003684_Cterm", "MF6_003686", "MF6_003843",
             "MF6_004284_Cterm", "MF6_004285", "MF6_004947", "MF6_006318"]
ABSENT = "#EDEDED"
NOTASSESSED = "#FFFFFF"   # white + x-hatch; never confuse with ABSENT grey
BRANCH = "#5A6068"
STAR = "#D55E00"
# a branch at or below this is treated as "no resolution" for the polytomy note
ZERO = 1e-7


def base_locus(locus):
    return locus[:-len("_Cterm")] if locus.endswith("_Cterm") else locus


def ring_label(locus):
    return (base_locus(locus) + " Cterm") if locus.endswith("_Cterm") else locus


# ------------------------------------------------------------------ tree source
# Prefer the IQ-TREE consensus (UFBoot, directly comparable to the chromosome-1
# tree's supports); fall back to FastTree, whose SH-like local supports are NOT
# bootstraps and must never be reported as though they were.
IQ, FT = BAC / "tables" / "iq_bac120.contree", BAC / "tables" / "tree_bac120.nwk"
if IQ.exists():
    TREEF, SUPPORT, SUPPORT_NOTE = IQ, "UFBoot", "1,000 ultrafast bootstrap replicates"
else:
    TREEF, SUPPORT, SUPPORT_NOTE = FT, "SH-like", ("FastTree SH-like local supports "
                                                   "— these are NOT bootstraps")
if not TREEF.exists():
    sys.exit(f"no bac120 tree found at {IQ} or {FT}")
print(f"tree: {TREEF.name}   support: {SUPPORT}")


# ------------------------------------------------------------------ data
def _acc(path):
    m = re.search(r"(GC[AF]_\d+\.\d+)", path or "")
    if m:
        return m.group(1)
    stem = Path(path or "").name
    for ext in (".fna", ".fa", ".fasta"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    return stem


# ANI: use the COMPLETE table (bac120/tables/mf6_ani_full.tsv), not
# rebuild/tables/MF6_ani_raw.tsv. The latter was computed against the
# DEREPLICATED set and covers only 553 of the 791 genomes here; 53 unmeasured
# genomes fall inside the MRCA of MF6's >=94% relatives, 37 of them B. sola.
# The full table is the SAME measurement -- same skani 0.3.2, same default flags,
# same whole-genome MF6 query -- re-run against the complete reference set. It
# reproduces the original on all 487 shared genomes to 2 d.p. (max diff 0.0000);
# that equality is asserted below rather than assumed.
ANIF = BAC / "tables" / "mf6_ani_full.tsv"
if not ANIF.exists():
    sys.exit(f"missing complete ANI table {ANIF} - run b9_ani_full.sh")
ani = {}
for r in csv.DictReader(open(ANIF), delimiter="\t"):
    if "MF6.fna" not in r["Query_file"]:
        continue
    other = _acc(r["Ref_file"])
    if other == "MF6":
        continue
    ani[other] = max(ani.get(other, 0.0), float(r["ANI"]))

_old = {}
for r in csv.DictReader(open(TAB / "MF6_ani_raw.tsv"), delimiter="\t"):
    a = _acc(r["Ref_file"])
    if a != "MF6":
        _old[a] = max(_old.get(a, 0.0), float(r["ANI"]))
_shared = set(_old) & set(ani)
_bad = [a for a in _shared if abs(ani[a] - _old[a]) >= 0.005]
assert not _bad, f"complete ANI table disagrees with the original on {len(_bad)} genomes: {_bad[:5]}"
print(f"ANI: {len(ani)} genomes measured (was {len(_old)}); "
      f"agrees with the original on all {len(_shared)} shared to 2 d.p.")

org, host = {}, {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    org[r["accession"]] = r["organism_name"]
    host[r["accession"]] = r["host_category"]
# the bac120 set brings in genomes the chr1 tables never saw; fill from its own
# metadata so a new tip is never silently drawn as "unknown"
_meta = json.load(open(BAC / "genome_metadata.json"))
for a, m in _meta.items():
    org.setdefault(a, m["organism"])

# 38 of this figure's 45 tips are absent from host_categories.tsv: they were
# removed by dereplication before Stage 10 ran, so the pipeline never classified
# them. Rather than invent a mapping or draw them all "unknown" -- 37 of them are
# B. sola from 'Maize rhizosphere', which is real, recorded provenance -- apply
# the SAME rules Stage 10 used (s10_hosts.py), in the same priority order:
# host_disease, then host, then isolation_source. Copied verbatim, not
# paraphrased; a re-implementation that drifted would silently recategorise.
NULL_TOKENS = {"", "na", "n/a", "none", "null", "unknown", "missing",
               "not collected", "not applicable", "not available",
               "not determined", "not provided", "unspecified", "-"}
HOST_RULES = [
    ("human_clinical", r"\bhomo sapiens\b|\bhuman\b"),
    ("fungus", r"rhizopus|fungus|fungal|mycelium|hypha|aspergillus|mucor"),
    ("amoeba_protist", r"dictyostelium|acanthamoeba|amoeba|paramecium"),
    ("insect", r"physopelta|insect|drosophila|beetle|aphid|termite|larva|hemiptera|coleoptera|bug\b"),
    ("animal", r"equus|horse|goat|sheep|ovis|capra|bovine|cattle|\bcow\b|\bpig\b|dog|canis|"
               r"felis|bird|chicken|gallus|mouse|mus musculus|\brat\b|rattus|monkey|reptile"),
    ("plant", r"oryza|\brice\b|sugarcane|arabidopsis|onion|zea mays|\bmaize\b|\bcorn\b|"
              r"glycine max|soybean|triticum|wheat|solanum|tomato|potato|vitis|musa|banana|"
              r"phaseolus|bean|medicago|mimosa|populus|eucalyptus|pine|pinus|bamboo|plant|tree"),
    ("soil", r"^soil$|^sand$"),
    ("water", r"^water$"),
    ("environmental_other", r"^environment$"),
]
SOURCE_RULES = [
    ("human_clinical", r"^human$|human |\bblood\b|sputum|\burine\b|cystic fibrosis|\bcf\b|"
                       r"clinical|patient|catheter|wound|abscess|nasal|bronch|tracheal|"
                       r"respiratory|melioidosis|sepsis|cepacia syndrome|swab|pus|cerebrospinal|"
                       r"bloodstream|blood stream|throat|lung|hospital|infection|surgical|milk"),
    ("rhizosphere", r"rhizosphere|rhizoplane|root nodule|nodule|rhizospheric"),
    ("plant", r"soybean|\bleaf\b|leaves|\broot\b|stem|seed|grain|panicle|\bplant\b|rice|"
              r"maize|corn|onion|sugarcane|phyllosphere|endophyt|fruit|flower|bark|tuber|"
              r"vegetable|crop|paddy|orchard|diseased"),
    ("fungus", r"fungus|fungal|hypha|mycelium|mushroom|lichen"),
    ("insect", r"insect|beetle|larva|gut of|termite|aphid"),
    ("animal", r"\bhorse\b|equine|goat|sheep|bovine|cattle|swine|veterinar|animal|glanders"),
    ("soil", r"\bsoil\b|sediment|compost|mud|\bfield\b|farmland|\bdust\b|sand"),
    ("water", r"\bwater\b|\briver\b|\blake\b|\bpond\b|marine|seawater|groundwater|aquatic|"
              r"stream|drinking|effluent|wastewater|sewage"),
    ("industrial", r"industrial|pharmaceutic|disinfectant|detergent|manufactur|antiseptic|"
                   r"mouthwash|cosmetic|food|dairy|beverage|fermentation|sludge"),
    ("environmental_other", r"environment|aerosol|\bair\b|atmospher|biofilm|bioreactor|enrich"),
]


def _clean(x):
    x = (x or "").strip().lower()
    return "" if x in NULL_TOKENS else x


def _match(rules, text):
    for cat, pat in rules:
        if re.search(pat, text):
            return cat
    return None


_derived = 0
for a, m in _meta.items():
    if a in host:
        continue
    h, src = _clean(m.get("host")), _clean(m.get("isolation_source"))
    cat = (_match(HOST_RULES, h) if h else None) or (_match(SOURCE_RULES, src) if src else None)
    host[a] = cat or "unknown"
    if cat:
        _derived += 1
print(f"host category derived by the Stage-10 rules for {_derived} genomes "
      f"the pipeline never classified")

c3 = {r["accession"]: (r["c3_present"] == "True")
      for r in csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t")}
# The 142 genomes that dereplication had removed were never typed. b12 typed them
# by s15b's method -- c3 := a secondary contig >=300 kb whose best skani hit
# against the KNOWN c3 contigs clears ANI >=95% and aligned fraction >=60% --
# which identifies c3 by similarity to already-typed replicons and so leaves the
# frozen classifier untouched. setdefault, never overwrite: a published call
# always wins over a newly computed one.
for r in csv.DictReader(open(BAC / "tables" / "c3_calls_new.tsv"), delimiter="\t"):
    c3.setdefault(r["accession"], r["c3_present"] == "True")
c3["MF6"] = True
# ASSESSED is load-bearing. 38 of this figure's 45 tips were removed by
# dereplication before the pC3-typing and toxin-search stages ran, so they have
# NO call at all. Drawing a missing call as "absent" -- which a bare dict .get()
# does -- would assert a measurement that was never made, and would do it exactly
# where the figure is making its point. Absent and never-assessed get different
# cells and different legend entries.
C3_ASSESSED = set(c3)

# Merged carrier table: the 773 retained calls verbatim, plus the 142 genomes
# dereplication had removed, searched with the same queries, the same thresholds
# and the same tier-1 rule (pident >=60, qcovs >=70).
# MF6_003686 IS A SPECIAL CASE and merging it from blastp would be wrong: it is a
# 120 aa query whose blastp coverage against the pyrodigal proteins caps at 65,
# below the tier-1 bar, because pyrodigal's ORF call for it is truncated. Run
# blastp alone and MF6 fails to carry its own locus. The retained pipeline used a
# dedicated tblastn against nucleotide, and so does the merge (b13).
carriers = defaultdict(set)
_assessed = defaultdict(set)
for r in csv.DictReader(open(BAC / "tables" / "toxin12_carriers_merged.tsv"),
                        delimiter="\t"):
    _assessed[r["locus"]].add(r["accession"])
    if r["tier1_carrier"] == "True":
        carriers[r["locus"]].add(r["accession"])

ORDER = sorted(RING_LOCI, key=lambda k: len(carriers[k]))
prot_col = {k: PROT_COLOR[i % len(PROT_COLOR)] for i, k in enumerate(ORDER)}

# ------------------------------------------------------------------ prune
tree = Phylo.read(str(TREEF), "newick")
tree.root_at_midpoint()
tree.ladderize()
alltips = {t.name: t for t in tree.get_terminals()}
N_INFERRED = len(alltips)          # before any display pruning
print(f"tips on the bac120 tree: {N_INFERRED}")
# strain names for the genomes whose accession alone is uninformative (the
# parent of the dropped mutant block is one of MF6's closest relatives, and
# "GCF_905400185.1" tells a reader nothing that "R-12632" does not tell better)
STRAIN = {r["parent_accession"]: r["parent_strain"] for r in td.lab_derivatives().values()}

if "MF6" not in alltips:
    sys.exit("MF6 is not on the bac120 tree")
# root -> MF6, so anc[-1-k] is the node k splits above the MF6 tip
anc = [tree.root] + list(tree.get_path(alltips["MF6"]))
if SPLITS_UP >= len(anc):
    sys.exit(f"MF6 has only {len(anc)-1} ancestors; --splits {SPLITS_UP} is above the root")
print("ancestors of MF6 and their clade sizes: "
      + ", ".join(f"{k}up={len(anc[-1-k].get_terminals())}"
                  for k in range(1, min(9, len(anc)))))
target = [t.name for t in anc[-1 - SPLITS_UP].get_terminals()]
print(f"clade {SPLITS_UP} splits above MF6: {len(target)} tips")

# THE SPLIT COUNT IS TAKEN ON THE FULL TREE, then the derivatives are pruned --
# never the other way round. Pruning first would collapse the single-child nodes
# it creates and silently renumber MF6's ancestors, so "6 splits up" would point
# at a different clade after the prune than before it.
dropped = td.drop_lab_derivatives(tree)
alltips = {t.name: t for t in tree.get_terminals()}
names = [n for n in target if n in alltips]
dropped_here = [n for n in target if n not in alltips]
clade = tree.common_ancestor([alltips[n] for n in names])
tips = clade.get_terminals()
assert {t.name for t in tips} == set(names), (
    "the clade is no longer exactly the retained target set after pruning")
print(f"drawing {len(names)} tips ({len(dropped_here)} isogenic lab derivatives dropped)")

# The guard that caught a real bug last time: a bad join made the "clade" come
# back as the whole tree. A close-relatives panel is tens of tips, not hundreds.
assert len(names) < 150, f"pruned clade has {len(names)} tips - selection is wrong"
assert "MF7" in names, "MF7 is not in the pruned clade - the point of the rebuild"

# What the topological cut costs relative to the old ANI cut, stated rather than
# hidden: tips this clade adds below the old 95% line, and >=95% relatives that
# fall OUTSIDE it (there must be none, or MF6's own neighbourhood is not a clade).
_in = set(names)
below = sorted((n for n in names if ani.get(n, 100.0) < ANI_MIN),
               key=lambda n: ani.get(n, 0.0))
outside = sorted((a for a, v in ani.items() if v >= ANI_MIN and a in alltips
                  and a not in _in), key=lambda a: -ani[a])
_anis = [ani[n] for n in names if n in ani]
print(f"ANI to MF6 across the clade: {min(_anis):.2f}-{max(_anis):.2f}%  "
      f"({len(below)} tips below {ANI_MIN:g}%)")
assert not outside, f"{len(outside)} genomes >={ANI_MIN}% ANI fall outside the clade: {outside[:5]}"

# Two different reasons a >=95% relative is not drawn, kept apart on purpose:
# never inferred at all, vs inferred and then dropped as an isogenic derivative.
_deriv = set(dropped)
missing = sorted((a for a, v in ani.items()
                  if v >= ANI_MIN and a not in alltips and a not in _deriv),
                 key=lambda a: -ani[a])
_deriv_hi = [a for a in _deriv if ani.get(a, 0.0) >= ANI_MIN]
print(f"relatives >={ANI_MIN}% not drawn: {len(missing)} never on the tree, "
      f"{len(_deriv_hi)} dropped as lab derivatives")

sup = [c.confidence for c in clade.get_nonterminals() if c.confidence is not None]

# ---- how much of this clade the markers actually resolve ---------------------
# Driven by the ALIGNMENT, not by branch lengths. Reading it off the tree was
# wrong, and wrong in a self-cancelling way: FastTree writes exact 0.0 for
# identical sequences but IQ-TREE clamps every branch to a floor of 2.55e-06, so
# a "<= 1e-7" test reported 43/45 unresolved tips on the FastTree tree and 0/45
# on the IQ-TREE one. The figure would have silently dropped its own resolution
# caveat at the moment it switched to the better tree. Sequence identity is a
# property of the data and does not move between tree methods.
_md5 = {r["accession"]: r["msa_md5"] for r in
        csv.DictReader(open(BAC / "tables" / "msa_seq_md5.tsv"), delimiter="\t")}
_by_hash = defaultdict(list)
for _t in tips:
    if _t.name in _md5:
        _by_hash[_md5[_t.name]].append(_t.name)
zero_tips = [n for v in _by_hash.values() if len(v) > 1 for n in v]
mf6_group = sorted(next((v for v in _by_hash.values() if "MF6" in v), []))
_ndist = len(_by_hash)
_floor = min((t.branch_length or 0.0) for t in tips)
_atfloor = sum(1 for t in tips if (t.branch_length or 0.0) <= _floor * 1.001)
_gap = [(t.name, k) for t in tips for k in RING_LOCI if t.name not in _assessed[k]]
_gapc3 = [t.name for t in tips if t.name not in C3_ASSESSED]
print(f"matrix coverage: {len(tips)*len(RING_LOCI) - len(_gap)}/{len(tips)*len(RING_LOCI)} "
      f"toxin cells, {len(tips)-len(_gapc3)}/{len(tips)} pC3 calls")
if _gap or _gapc3:
    print(f"  STILL UNASSESSED: {len(_gap)} toxin cells, {len(_gapc3)} pC3 calls")
print(f"distinct marker sequences among the {len(tips)} tips: {_ndist}")
print(f"tips sharing a sequence with another: {len(zero_tips)}/{len(tips)}")
print(f"terminal branches at the tree floor ({_floor:.3g}): {_atfloor}/{len(tips)}")
print(f"MF6's identical group: {mf6_group}")

# ------------------------------------------------------------------ layout
xpos, ypos = {}, {}
xpos[id(clade)] = 0.0
for cl in clade.find_clades(order="preorder"):
    for c in cl.clades:
        xpos[id(c)] = xpos[id(cl)] + max(c.branch_length or 0.0, 0.0)
for i, t in enumerate(tips):
    ypos[id(t)] = float(i)
for cl in clade.find_clades(order="postorder"):
    if cl.clades:
        ypos[id(cl)] = float(np.mean([ypos[id(c)] for c in cl.clades]))
XMAX = max(xpos[id(t)] for t in tips)

# Height follows the tip count so the row pitch stays legible: the header and
# footer blocks are fixed in INCHES (they hold a fixed number of text lines), so
# every figure-fraction y below is computed from them rather than hard-coded --
# hard-coded fractions silently slide the captions into the panel as soon as the
# clade grows. Row pitch is capped so an 84-tip clade does not produce a 12-inch
# figure, and the tip fonts follow the pitch.
HDR_IN, FTR_IN = 1.264, 1.443      # as in the 45-tip version (7.8 in tall)
HMAX = 11.0
PITCH = min(0.113, max(0.070, (HMAX - HDR_IN - FTR_IN) / len(tips)))
FIGH = HDR_IN + FTR_IN + PITCH * len(tips)
TIPFS = float(np.clip(PITCH * 72 * 0.78, 4.6, 6.3))
ANIFS = TIPFS - 0.5
# headroom above the first tip for the rotated matrix column labels: a fixed
# 0.75 in of paper, expressed in rows, so it does not shrink with the pitch
HEAD_ROWS = max(6.4, 0.75 / PITCH)


def _ytop(inches):      # figure fraction, `inches` down from the top edge
    return 1.0 - inches / FIGH


def _ybot(inches):      # figure fraction, `inches` up from the bottom edge
    return inches / FIGH


print(f"figure: {FIGH:.2f} in tall, {PITCH:.3f} in/tip, tip labels {TIPFS:.1f} pt")
fig = plt.figure(figsize=(9.8, FIGH))
gsp = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.0],
                       left=0.035, right=0.985,
                       top=_ytop(HDR_IN), bottom=_ybot(FTR_IN),
                       wspace=0.02)
axT = fig.add_subplot(gsp[0, 0])
axM = fig.add_subplot(gsp[0, 1], sharey=axT)

segs = []
for cl in clade.find_clades(order="preorder"):
    if not cl.clades:
        continue
    ys = [ypos[id(c)] for c in cl.clades]
    segs.append([(xpos[id(cl)], min(ys)), (xpos[id(cl)], max(ys))])
    for c in cl.clades:
        segs.append([(xpos[id(cl)], ypos[id(c)]), (xpos[id(c)], ypos[id(c)])])
axT.add_collection(LineCollection(segs, colors=BRANCH, linewidths=0.9, zorder=2))

# Print a support only if it means something. IQ-TREE fully resolves the
# identical-sequence polytomies in this clade, and its own bootstrap says those
# splits are worthless: across the whole tree the 453 splits FastTree does NOT
# recover have median UFBoot 19, against median 95 for the 335 both methods find.
# Printing 40 two-digit numbers in the teens would decorate the panel with noise
# and imply a resolution the data does not support.
SUP_MIN = 50
_printed = 0
for cl in clade.get_nonterminals():
    if cl.confidence is not None and SUP_MIN <= cl.confidence < 100 and cl is not clade:
        axT.text(xpos[id(cl)] - XMAX * 0.006, ypos[id(cl)] + 0.30,
                 f"{cl.confidence:.0f}", fontsize=5.4, ha="right", va="center",
                 color="#777777", zorder=5)
        _printed += 1
_low = sum(1 for cl in clade.get_nonterminals()
           if cl.confidence is not None and cl.confidence < SUP_MIN)
print(f"support labels printed: {_printed}   suppressed as <{SUP_MIN}: {_low}")

LX = XMAX * 1.035
for t in tips:
    y = ypos[id(t)]
    inmf6 = t.name in mf6_group
    axT.plot([xpos[id(t)], LX - XMAX * 0.008], [y, y], lw=0.4, ls=(0, (1, 2)),
             color="#C8C8C8", zorder=1)
    sp = org.get(t.name, t.name)
    sp = sp.replace("Burkholderia ", "B. ").replace("Paraburkholderia ", "P. ")
    if t.name in ("MF6", "MF7"):
        lab = f"B. sola {t.name}"
    elif t.name in STRAIN:
        lab = f"{sp} {STRAIN[t.name]}   {t.name}"
    else:
        lab = f"{sp}   {t.name}"
    star = t.name in ("MF6", "MF7")
    axT.text(LX, y, lab, fontsize=TIPFS, va="center", ha="left",
             color=STAR if star else "#222222",
             fontweight="bold" if star else "normal", zorder=5)
    if star:
        # BOTH filled now: on this tree MF7 is inferred, not grafted.
        axT.scatter([LX - XMAX * 0.030], [y], s=52, marker="*", c=STAR,
                    edgecolors="black", linewidths=0.4, zorder=6)
    a = ani.get(t.name)
    if a is not None:
        axT.text(LX + XMAX * 0.86, y, f"{a:.2f}", fontsize=ANIFS, va="center",
                 ha="right", color="#666666", zorder=5)

# bracket the unresolved group so the polytomy is explicit rather than implied
if mf6_group:
    # PER-TIP MARKS, not a bracket. A bracket spans a contiguous y-range, and on
    # the IQ-TREE topology this group is NOT contiguous -- GCF_018417615.1 sits
    # between its members because IQ-TREE arbitrarily resolved the polytomy. A
    # bracket would have asserted that LAS2 is one of the identical genomes.
    ys = [ypos[id(t)] for t in tips if t.name in mf6_group]
    bx = XMAX * 1.945
    for y in ys:
        axT.plot([bx - XMAX * 0.014, bx + XMAX * 0.014], [y, y], lw=2.0,
                 color=STAR, solid_capstyle="butt", zorder=6)
    axT.text(bx + XMAX * 0.030, float(np.mean(ys)),
             f"these {len(mf6_group)} are identical\nacross all 5,010 columns",
             fontsize=5.3, va="center", ha="left", color=STAR, rotation=90,
             linespacing=1.25, zorder=6)

axT.text(LX + XMAX * 0.86, -HEAD_ROWS * 0.18, "ANI to\nMF6 (%)", fontsize=5.8, ha="right",
         va="center", color="#666666", linespacing=1.3)
axT.set_xlim(-XMAX * 0.02, XMAX * 2.10)
axT.set_ylim(len(tips) - 0.4, -HEAD_ROWS)
axT.axis("off")

_n = 10 ** np.floor(np.log10(XMAX / 3))
step = _n * max(1, round((XMAX / 3) / _n))
axT.plot([0, step], [len(tips) - 0.75] * 2, lw=1.2, color="#000000")
axT.text(step / 2, len(tips) - 1.05, f"{step:g} subs/site", ha="center",
         va="bottom", fontsize=6)

# ------------------------------------------------------------------ matrix
COLS = ([("pC3", None, figstyle.C3_COLORS[True]),
         ("isolation source", None, "#BBBBBB")]
        + [(ring_label(k), k, prot_col[k]) for k in ORDER])
for j, (lab, key, col) in enumerate(COLS):
    for t in tips:
        y = ypos[id(t)]
        hatch = None
        if key is None and lab == "pC3":
            if t.name not in C3_ASSESSED:
                fc, hatch = NOTASSESSED, "xxx"
            else:
                fc = figstyle.C3_COLORS[True] if c3[t.name] else ABSENT
        elif key is None:
            h = host.get(t.name, "unknown")
            fc = figstyle.HOST_COLORS.get(h, ABSENT)
            if h == "unknown":
                hatch = "///"
        elif t.name not in _assessed[key]:
            fc, hatch = NOTASSESSED, "xxx"
        else:
            fc = col if t.name in carriers[key] else ABSENT
        axM.add_patch(Rectangle((j + 0.06, y - 0.42), 0.88, 0.84,
                                facecolor=fc, hatch=hatch,
                                edgecolor="#B0B0B0" if hatch else "none",
                                linewidth=0.3))
    axM.text(j + 0.5, -0.75, lab, fontsize=5.9, rotation=90, ha="center",
             va="bottom", color="#222222")
axM.plot([2.0, 2.0], [-0.55, len(tips) - 0.45], color="#FFFFFF", lw=2.4, zorder=4)
axM.set_xlim(0, len(COLS))
axM.axis("off")

# ------------------------------------------------------------------ labels
# Every number quoted below is over the ASSESSED subset, with the denominator
# stated. "7/45 carry pC3" would be false -- 38 of the 45 were never tested.
_ass_c3 = [t for t in tips if t.name in C3_ASSESSED]
_ass_tox = [t for t in tips if all(t.name in _assessed[k] for k in RING_LOCI)]
n_pc3 = sum(1 for t in _ass_c3 if c3[t.name])
n_loci = {t.name: sum(1 for k in RING_LOCI if t.name in carriers[k]) for t in _ass_tox}
n_unassessed = len(tips) - len(_ass_c3)
_full = [t.name for t in tips if n_loci.get(t.name, 0) == len(RING_LOCI)]
_partial = sorted(((n_loci[t.name], t.name) for t in tips
                   if n_loci.get(t.name, 0) < len(RING_LOCI)))
print(f"carrying all {len(RING_LOCI)} loci: {len(_full)}/{len(tips)}")
print("  not complete: " + ", ".join(f"{n}/10 {a}" for n, a in _partial))

# "hidden by dereplication" is COUNTED, not assumed: a genome is retained iff it
# carried a published pC3 call before this rebuild. Hard-coding the number (it
# was 3 when the panel held 45 tips) silently lies as soon as the clade changes.
RETAINED = {r["accession"] for r in
            csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t")}
n_ret = sum(1 for t in tips if t.name in RETAINED)
n_hidden_full = sum(1 for a in _full if a not in RETAINED)
print(f"tips retained through dereplication: {n_ret}/{len(tips)}; "
      f"complete-set carriers hidden by it: {n_hidden_full}/{len(_full)}")

def _short(a):
    if a in STRAIN:
        return STRAIN[a]
    nm = (org.get(a, a) or a).replace("Burkholderia ", "B. ").split()
    return nm[-1] if nm and nm[-1] not in ("sola", "sp.", "cenocepacia") else a


_others = sorted((a for a in _full if a not in ("MF6", "MF7")),
                 key=lambda a: -ani.get(a, 0.0))
fig.text(0.035, _ytop(0.14),
         f"The complete ten-locus set is not unique to MF6 and MF7 — {len(_others)} other\n"
         f"genome{'s' if len(_others) != 1 else ''} of the {len(tips)} here carr{'y' if len(_others) != 1 else 'ies'} all ten: "
         + ", ".join(f"{_short(a)}" for a in _others),
         fontsize=10.5, fontweight="bold", va="top", linespacing=1.35)

_res = (f"{len(mf6_group)} genomes ({', '.join(mf6_group)}) are byte-identical across all "
        f"5,010 masked columns" if mf6_group else "no identical group found")
_spnames = sorted({" ".join((org.get(t.name, "") or "").split()[:2]) for t in tips} - {""})
_sp = ", ".join(n.replace("Burkholderia ", "B. ") for n in _spnames)
fig.text(0.035, _ytop(0.74),
         f"Pruned from the {N_INFERRED}-tip bac120 marker tree (GTDB-Tk r232, 120 single-copy markers, 5,010 masked AA columns; not dereplicated).\n"
         f"Tips are EVERY descendant of the node {SPLITS_UP} splits above MF6 — a clade, not an ANI cut — spanning {min(_anis):.2f}–{max(_anis):.2f}% ANI to MF6 "
         f"({len(below)} of\n{len(tips)} fall below the {ANI_MIN:g}% line the previous version used): {_sp}. "
         f"{len(dropped_here)} tips of that clade are NOT drawn: they are isogenic transposon mutants of\n"
         f"B. sola R-12632 (PRJEB40633; Depoorter et al. 2021 AEM 87:e01169-21), one strain sequenced 38 times, and the parent is kept in their place. "
         f"Support is {SUPPORT} ({SUPPORT_NOTE}); {_low} of {_low + _printed + 1} internal\nnodes "
         f"score below 50 and their values are suppressed rather than printed as if meaningful. "
         f"pC3 is present in {n_pc3}/{len(_ass_c3)} tips. Every cell is measured: the "
         f"{len(tips) - n_ret} genomes that\ndereplication had removed were newly typed for pC3 "
         f"(skani vs known c3 replicons) and searched for all ten loci with the retained "
         f"queries and thresholds.",
         fontsize=6.3, color="#444444", va="top", linespacing=1.5)

fig.text(0.035, _ybot(FTR_IN - 0.03),
         f"RESOLUTION LIMIT — read this figure alongside Figure 12, not instead of it: {_res}.\n"
         f"only {_ndist} of the {len(tips)} tips have a distinct marker sequence: 120 genes at 5,010 positions cannot separate genomes >99.5% identical, whereas the\n"
         f"1,337,868 bp chromosome-1 alignment behind Figure 12 does. What the rebuild buys: MF7 is INFERRED here, not grafted (all 120 markers recovered\n"
         f"from its 63-contig draft, which the chromosome-1 alignment dropped at 93.96% gaps), and MF6's closest relatives are ordinary tips rather than a\n"
         f"footnote. The matrix is now complete — {len(tips)*len(RING_LOCI)}/{len(tips)*len(RING_LOCI)} toxin cells and {len(tips)}/{len(tips)} pC3 calls measured, not inferred. MF6_003686 needed a dedicated tblastn:\n"
         f"its 120 aa blastp coverage caps at 65 against the pyrodigal ORFs, below the tier-1 bar, so blastp alone scores MF6 as not carrying its own locus.",
         fontsize=5.9, color=STAR, va="top", linespacing=1.5)

hosts_here = [h for h in figstyle.HOST_COLORS
              if any(host.get(t.name) == h for t in tips)]
leg = fig.legend(handles=[Patch(facecolor=figstyle.C3_COLORS[True], label="pC3 present"),
                          Patch(facecolor=ABSENT, label="absent (measured)")]
                 + ([Patch(facecolor=NOTASSESSED, hatch="xxx", edgecolor="#999999",
                           linewidth=0.3, label="never assessed")]
                    if (_gap or _gapc3) else [])
                 + [Patch(facecolor=figstyle.HOST_COLORS[h],
                          edgecolor="#B0B0B0" if h == "unknown" else "#999999",
                          hatch="///" if h == "unknown" else None, linewidth=0.3,
                          label="source not recorded" if h == "unknown"
                          else h.replace("_", " ")) for h in hosts_here],
                 loc="lower center", ncol=8, frameon=False, fontsize=6.2,
                 handlelength=1.2, columnspacing=1.1, bbox_to_anchor=(0.5, _ybot(0.016)),
                 title="isolation source · toxin columns use the Figure 11 locus colours",
                 title_fontsize=6.8)
leg.get_title().set_fontweight("bold")

ok = figstyle.save(fig, str(FIG / "fig12b_close_relatives_bac120"))
print("OK" if ok else "FAILED vector-text check")
