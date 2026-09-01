#!/usr/bin/env python3
"""Figure 13 - genus-wide bac120 single-copy marker phylogeny of Burkholderia s.l.

The taxonomic backbone Moshe asked for: every available Burkholderia sensu lato
genome, NOT dereplicated, with pC3 presence and isolation source mapped onto it.

HOW IT DIFFERS FROM FIGURE 4, and why both exist:
  Figure 4  chromosome-1 core genes, 763 dereplicated genomes, 1,337,868 bp.
            High resolution; one representative per near-clone cluster.
  Figure 13 120 GTDB single-copy markers, 791 genomes inferred / 754 drawn,
            5,010 masked AA columns.
            Every genome present; far LOWER resolution -- see the caveat below.

THE RESOLUTION CAVEAT IS THE HEADLINE, not a footnote. Measured directly on the
alignment: only 448 of the 791 inferred sequences are DISTINCT. 427 genomes (54.0%) are
byte-identical to at least one other across all 5,010 columns, in 84 groups, the
largest of 70 genomes. IQ-TREE resolves those polytomies anyway -- its branch
floor is 2.55e-06, not 0 -- but its own bootstrap disowns the result: of the 788
internal splits, the 453 that FastTree does NOT recover have median UFBoot 19,
against median 95 for the 335 both methods find.
So: this tree is the right instrument for "which LINEAGES carry pC3, and from
what sources", and the wrong one for "how are these near-clones related". For
the latter use Figure 4 or ANI.

Every tip carries measured data. The 142 genomes dereplication had removed were
newly typed for pC3 (b12, skani vs known c3 replicons -- the frozen classifier is
untouched) and classified for isolation source with Stage 10's own rules.

Rings, inner to outer:  1. pC3 presence   2. isolation source
"""
import collections
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import radialtree as rt  # noqa: E402
import tree_display as td  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
BAC = D.parent / "bac120"
FIG.mkdir(exist_ok=True)

IQ, FT = BAC / "tables" / "iq_bac120.contree", BAC / "tables" / "tree_bac120.nwk"
if IQ.exists():
    TREE, SUPPORT, SUPPORT_NOTE = IQ, "UFBoot", "1,000 ultrafast bootstrap replicates"
else:
    TREE, SUPPORT, SUPPORT_NOTE = FT, "SH-like", "FastTree SH-like supports, NOT bootstraps"
if not TREE.exists():
    sys.exit(f"no bac120 tree at {IQ} or {FT}")
print(f"tree: {TREE.name}   support: {SUPPORT}")

tree = Phylo.read(str(TREE), "newick")
tree.root_at_midpoint()
tree.ladderize()
N_INFERRED = len(tree.get_terminals())     # before any display pruning
# Drop isogenic laboratory derivatives (see tree_display.drop_lab_derivatives):
# BioProject PRJEB40633 put one B. sola strain, R-12632, on this tree 38 times --
# the parent plus 37 transposon mutants of it, each inheriting the parent's
# "Maize rhizosphere, Italy, 1996" metadata. Drawn, they are 4.7% of the rim and
# a 34-genome identical-sequence block, all pC3-positive, sitting on MF6's
# doorstep: they would inflate both the pC3 fraction and the largest-identical-
# group statistic with a single genome. The parent stays; the tree and the
# alignment are untouched, the mutants are simply not drawn.
DROPPED = td.drop_lab_derivatives(tree)
tips_all = tree.get_terminals()
print(f"tips: {len(tips_all)} drawn ({N_INFERRED} inferred, "
      f"{len(DROPPED)} isogenic lab derivatives dropped)")

# ---------------- trait tables ----------------
# organism / host, with the same Stage-10 fallback fig12b uses for the genomes
# the pipeline never classified. Copied rules, not paraphrased ones.
host, org = {}, {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    host[r["accession"]] = r["host_category"]
    org[r["accession"]] = r["organism_name"]
_meta = json.load(open(BAC / "genome_metadata.json"))
for a, m in _meta.items():
    org.setdefault(a, m["organism"])

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
    _derived += bool(cat)
print(f"isolation source derived by the Stage-10 rules for {_derived} genomes")

# pC3: retained calls verbatim, then the newly typed ones. setdefault, never
# overwrite -- a published call always wins.
c3_present = {}
for r in csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t"):
    c3_present[r["accession"]] = (r["c3_present"] == "True")
_new = 0
for r in csv.DictReader(open(BAC / "tables" / "c3_calls_new.tsv"), delimiter="\t"):
    if r["accession"] not in c3_present:
        c3_present[r["accession"]] = (r["c3_present"] == "True"); _new += 1
c3_present["MF6"] = True
print(f"pC3 calls: {len(c3_present)} total, {_new} newly measured for this figure")

_unc3 = [t.name for t in tips_all if t.name not in c3_present]
if _unc3:
    sys.exit(f"{len(_unc3)} tips have no pC3 call, e.g. {_unc3[:3]} - "
             "a missing call must never be drawn as an absence")

# ---------------- how much the markers resolve ----------------
# From the ALIGNMENT, not from branch lengths: IQ-TREE clamps branches to a floor
# of 2.55e-06 rather than 0, so a branch-length test silently reports full
# resolution on this tree while reporting none on the FastTree one.
_md5 = {r["accession"]: r["msa_md5"] for r in
        csv.DictReader(open(BAC / "tables" / "msa_seq_md5.tsv"), delimiter="\t")}
_grp = collections.defaultdict(list)
for t in tips_all:
    if t.name in _md5:
        _grp[_md5[t.name]].append(t.name)
N_DISTINCT = len(_grp)
N_DUP = sum(len(v) for v in _grp.values() if len(v) > 1)
print(f"distinct marker sequences: {N_DISTINCT}/{len(tips_all)}   "
      f"genomes sharing one: {N_DUP} ({100*N_DUP/len(tips_all):.1f}%)   "
      f"largest identical group: {max(len(v) for v in _grp.values())}")

# ---------------- ingroup / outgroup / near-clone wedges ----------------
print("display conventions:")
ingroup, outgroup = td.ingroup_root(tree, c3_present, org)
out_tips = outgroup.get_terminals()
nearclones = td.nearclone_clades(tree, c3_present, org, skip=[outgroup])
_hidden = {t.name for cl in nearclones for t in cl.get_terminals()}
tips = [t for t in tree.get_terminals()
        if t.name not in {x.name for x in out_tips} and t.name not in _hidden]

n_pres = sum(1 for t in tips_all if c3_present.get(t.name))
print(f"pC3 present on tree: {n_pres}/{len(tips_all)}")

# ---------------- draw ----------------
fig = plt.figure(figsize=(8.4, 9.1))
ax = fig.add_axes([0.02, 0.190, 0.96, 0.735])
ax.set_aspect("equal"); ax.axis("off")

OG_WEIGHT, NC_WEIGHT = 22.0, 10.0
_collapse = {outgroup: OG_WEIGHT}
_collapse.update({cl: NC_WEIGHT for cl in nearclones})
lay = rt.RadialLayout(tree, start_deg=92, extent_deg=352, collapse=_collapse)
_slots = len(tips) + OG_WEIGHT + NC_WEIGHT * len(nearclones)
print(f"drawn elements: {len(tips)} tips + {len(nearclones)} near-clone wedges "
      f"+ 1 outgroup wedge (from {len(tips_all)} tips)")
print(f"angular pitch: {352.0/_slots:.4f} deg/slot "
      f"({352.0/_slots/(352.0/len(tips_all)):.2f}x the uncollapsed pitch)")

state = {}
for t in tree.get_terminals():
    state[id(t)] = host.get(t.name, "unknown")
for cl in tree.get_nonterminals(order="postorder"):
    s = {state[id(c)] for c in cl.clades}
    state[id(cl)] = s.pop() if len(s) == 1 else None

GREY = "#C8C8C8"


def colour_of(node):
    st = state.get(id(node))
    return figstyle.HOST_COLORS.get(st, GREY) if st else GREY


segs, cols = [], []
for cl in tree.get_nonterminals(order="preorder"):
    if id(cl) not in lay.radius:
        continue
    ra = lay.radius[id(cl)]
    angs = [lay.angle[id(c)] for c in cl.clades if id(c) in lay.radius]
    if len(angs) > 1:
        xs, ys = rt.arc_points(ra, min(angs), max(angs),
                               n=max(4, int(np.rad2deg(max(angs) - min(angs)))))
        segs.append(list(zip(xs, ys))); cols.append(colour_of(cl))
    for c in cl.clades:
        if id(c) not in lay.radius:
            continue
        rc, ac = lay.radius[id(c)], lay.angle[id(c)]
        segs.append([rt.polar_to_xy(ra, ac), rt.polar_to_xy(rc, ac)])
        cols.append(colour_of(c))
ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.55, zorder=2))
n_painted = sum(1 for c in cols if c != GREY)
print(f"branch segments painted by source: {n_painted}/{len(cols)} "
      f"({100*n_painted/len(cols):.0f}%)")

R0 = 1.04
_ext = [[rt.polar_to_xy(lay.tip_radius(t), lay.tip_angle(t)),
         rt.polar_to_xy(R0, lay.tip_angle(t))] for t in tips]
ax.add_collection(LineCollection(_ext, colors="#BBBBBB", linewidths=0.22,
                                 linestyles=(0, (1, 2)), zorder=1))

rt.draw_ring(ax, lay, {t.name: c3_present.get(t.name, False) for t in tips},
             figstyle.C3_COLORS, R0, R0 + 0.055, default="#D9D9D9")
rt.draw_ring(ax, lay, {t.name: host.get(t.name, "unknown") for t in tips},
             figstyle.HOST_COLORS, R0 + 0.068, R0 + 0.123, default="#EEEEEE")

OG_CAP = 0.86
for cl, a0, a1, r_stem, r_crown in lay.collapsed_wedges():
    truncated = r_crown > 1.0
    r_cap = OG_CAP if truncated else r_crown
    td.draw_wedge(ax, rt, a0, a1, r_stem, r_cap, break_marker=truncated)
    _la = np.linspace(a0, a1, 7)[1:-1]
    ax.add_collection(LineCollection(
        [[rt.polar_to_xy(r_cap, _a), rt.polar_to_xy(R0, _a)] for _a in _la],
        colors="#BBBBBB", linewidths=0.22, linestyles=(0, (1, 2)), zorder=1))
    td.ring_cell(ax, a0, a1, R0, R0 + 0.055, figstyle.C3_COLORS[False])
    _h = {host.get(t.name, "unknown") for t in cl.get_terminals()}
    _hc = figstyle.HOST_COLORS.get(next(iter(_h)), "#EEEEEE") if len(_h) == 1 else GREY
    td.ring_cell(ax, a0, a1, R0 + 0.068, R0 + 0.123, _hc)

    am = 0.5 * (a0 + a1)
    if cl is outgroup:
        _og = collections.Counter(
            (org.get(t.name, "?").split()[0]) for t in cl.get_terminals())
        _txt = " · ".join(f"$\\it{{{g}}}$ {n}" for g, n in _og.most_common(4))
        # +0.30, not +0.215: at the smaller radius this caption crossed the
        # rotated "B. pseudomallei (263)" wedge label next to it.
        lx, ly = rt.polar_to_xy(R0 + 0.300, am)
        ax.text(lx, ly, f"outgroup, {len(cl.get_terminals())} genomes\n{_txt}\n"
                        f"collapsed; true depth {r_crown:.1f}× this radius · no pC3",
                fontsize=5.9, color="#333333",
                ha="right" if np.cos(am) < 0 else "left",
                va="bottom" if np.sin(am) > 0 else "top", zorder=9,
                linespacing=1.35)
    else:
        deg = np.rad2deg(am) % 360.0
        flip = 90.0 < deg < 270.0
        lx, ly = rt.polar_to_xy(R0 + 0.135, am)
        ax.text(lx, ly, td.clade_label(cl, org, italic=True, abbrev=True),
                fontsize=6.0, color="#333333", zorder=9,
                rotation=deg + 180 if flip else deg,
                rotation_mode="anchor",
                ha="right" if flip else "left", va="center")

# MF6 and MF7 both get a star here. On the chromosome-1 tree MF7 could not be
# drawn at all -- its 63-contig draft was dropped from that alignment at 93.96%
# gaps -- and this is the figure where it is a real inferred tip: all 120 markers
# recovered, byte-identical marker counts to MF6.
_starred = 0
for t in tips:
    if t.name in ("MF6", "MF7"):
        a = lay.tip_angle(t)
        x, y = rt.polar_to_xy(R0 + 0.185, a)
        ax.scatter([x], [y], s=64, marker="*", c="#D55E00",
                   edgecolors="black", linewidths=0.5, zorder=10)
        xt, yt = rt.polar_to_xy(R0 + 0.235, a)
        ax.text(xt, yt, t.name, fontsize=7.5, fontweight="bold", color="#D55E00",
                ha="left" if np.cos(a) > 0 else "right", va="center", zorder=10)
        _starred += 1
if _starred < 2:
    print(f"NOTE: only {_starred} of MF6/MF7 drawn as tips "
          "(the other is inside a collapsed near-clone wedge)")

lim = R0 + 0.52
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

depth = lay.max_r
_n = 10 ** np.floor(np.log10(depth / 4))
step = _n * max(1, round((depth / 4) / _n))
ax.plot([-lim + 0.04, -lim + 0.04 + step / depth], [-lim + 0.10, -lim + 0.10],
        lw=1.2, color="#000000")
ax.text(-lim + 0.04 + step / (2 * depth), -lim + 0.14, f"{step:g} subs/site",
        ha="center", va="bottom", fontsize=7)

fig.text(0.5, 0.985, "pC3 across the bac120 single-copy marker phylogeny of "
         "$\\it{Burkholderia}$ sensu lato",
         ha="center", va="top", fontsize=10.5, fontweight="bold")
fig.text(0.5, 0.958,
         f"{len(tips_all)} genomes, NOT dereplicated — every available complete/chromosome-level "
         f"assembly with isolation metadata, less {len(DROPPED)} isogenic transposon\n"
         f"mutants of B. sola R-12632 (PRJEB40633) whose parent is kept. 120 GTDB markers, 5,010 masked AA\n"
         f"columns, IQ-TREE LG+F+R4 with {SUPPORT_NOTE}. pC3 present in {n_pres} "
         f"({100*n_pres/len(tips_all):.0f}%); the {_new} genomes dereplication had removed were "
         f"newly typed, none inherited.\n"
         f"Rooted on the ingroup stem; the {len(out_tips)}-genome outgroup and "
         f"{len(nearclones)} near-clone blocks ({len(_hidden)} genomes) are collapsed. "
         f"Branches painted by source where all descendants agree.",
         ha="center", va="top", fontsize=7.2, color="#444444")

fig.text(0.5, 0.052,
         f"RESOLUTION: only {N_DISTINCT} of the {len(tips_all)} drawn genomes have a DISTINCT marker "
         f"sequence — {N_DUP} ({100*N_DUP/len(tips_all):.0f}%) are byte-identical to another across all "
         f"5,010 columns, the largest group being {max(len(v) for v in _grp.values())} genomes.\n"
         f"IQ-TREE resolves those polytomies anyway (its branch floor is 2.55e-06, not 0) but its own "
         f"bootstrap disowns them: on the full {N_INFERRED}-tip inference the 453 splits FastTree does not recover have median UFBoot 19,\n"
         f"against 95 for the 335 both methods find. Read this tree for WHICH LINEAGES carry pC3 and "
         f"from what sources — not for relationships among near-clones, where Figure 4 or ANI applies.",
         ha="center", va="top", fontsize=5.9, color="#B35A00", linespacing=1.5)

h1 = [Patch(facecolor=figstyle.C3_COLORS[True], label="pC3 present"),
      Patch(facecolor=figstyle.C3_COLORS[False], label="pC3 absent")]
present_hosts = [h for h in figstyle.HOST_COLORS
                 if any(host.get(t.name) == h for t in tips_all)]
h2 = [Patch(facecolor=figstyle.HOST_COLORS[h], label=h.replace("_", " "))
      for h in present_hosts]
h4 = [Line2D([0], [0], color=GREY, lw=1.6, label="branch: sources disagree"),
      Patch(facecolor="#E8E4DC", edgecolor="#7A7A7A",
            label="collapsed clade (pC3-negative)")]
leg1 = fig.legend(handles=h1 + h2 + h4, loc="lower center", ncol=8,
                  frameon=False, fontsize=6.5, handlelength=1.2,
                  columnspacing=1.1, labelspacing=0.32,
                  bbox_to_anchor=(0.5, 0.088),
                  title="pC3 presence  ·  isolation source",
                  title_fontsize=7)
leg1.get_title().set_fontweight("bold")

ok = figstyle.save(fig, str(FIG / "fig13_genus_tree_bac120"))
print("OK" if ok else "FAILED vector-text check")
