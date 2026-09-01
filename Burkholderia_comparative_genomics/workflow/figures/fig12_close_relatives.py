#!/usr/bin/env python3
"""Figure 12 - MF6 and its closest relatives, pruned out of the chromosome-1 tree.

Michal asked for "a small tree of the very closely related strains only". This is
that tree, and it is a PRUNE of tables/chr1_core.treefile -- not a new inference.
Every branch length and every UFBoot value is carried over unchanged from the
304-tip analysis, so this figure and Figure 11 cannot disagree.

Which tips: MF6 plus every genome on the 304-tip tree at >=94% ANI to MF6
(tables/MF6_ani_raw.tsv). That set turns out to be EXACTLY monophyletic -- the
MRCA clade of those 18 tips contains those 18 tips and nothing else -- so no
intervening lineage is being hidden by the selection. That is worth stating,
because pruning to a hand-picked list of tips is otherwise a good way to make
unrelated strains look like each other's nearest neighbours. The same set comes
out at >=93%; at >=92% it grows to 33 tips (30 selected + 3 intervening) and at
>=91% it explodes to 101, so 94% is both the natural and the conservative cut.

CAVEAT, and it is the important one: the single closest relative of MF6,
GCF_016899425.1 (Burkholderia sp. MS389, 99.71% ANI), is NOT here. It was
removed by the dereplication step that produced the 304-tip tree, so it is not
available to prune to. GCF_053209605.1 (96.11%) is absent for the same reason.
Both are present in the toxin search and in Figure 11 panel B. Putting them on a
tree requires re-running the phylogeny from genomes, which this figure does not
do. The caption must say so.

Columns, left to right: pC3 presence, isolation source, then the ten candidate
toxin / immunity loci in the same order and the same colours as Figure 11's
rings, so the two figures read as one story. As in Figure 11, the two RHS loci
are represented by their C-terminal tip, not the full-length scaffold.

Inputs: tables/chr1_core.treefile, MF6_ani_raw.tsv, host_categories.tsv,
        c3_calls_all_genomes.tsv, toxin12_carriers.tsv, toxin12_locus_summary.tsv
"""
import csv
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
from Bio.Phylo.BaseTree import Clade  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

ANI_MIN = 94.0

# kept byte-identical to fig11_toxin_tree.py so a locus is the same colour in
# both figures; if you change one, change the other
PROT_COLOR = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
              "#DDCC77", "#CC6677", "#882255", "#AA4499", "#EE7733"]
RING_LOCI = ["MF6_001079", "MF6_002734", "MF6_003184",
             "MF6_003684_Cterm", "MF6_003686", "MF6_003843",
             "MF6_004284_Cterm", "MF6_004285", "MF6_004947", "MF6_006318"]
ABSENT = "#EDEDED"
BRANCH = "#5A6068"
STAR = "#D55E00"


def base_locus(locus):
    return locus[:-len("_Cterm")] if locus.endswith("_Cterm") else locus


def ring_label(locus):
    return (base_locus(locus) + " Cterm") if locus.endswith("_Cterm") else locus


# ------------------------------------------------------------------ data
def _acc(path):
    """File path -> accession. Falls back to the bare stem so in-house genomes
    (MF6, MF7) come back as "MF7", not "genomes/MF7.fna"; the original regex
    only knew NCBI accessions and leaked the path into tip lookups and captions."""
    m = re.search(r"(GC[AF]_\d+\.\d+)", path or "")
    if m:
        return m.group(1)
    stem = Path(path or "").name
    for ext in (".fna", ".fa", ".fasta"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
    return stem


ani = {}
for r in csv.DictReader(open(TAB / "MF6_ani_raw.tsv"), delimiter="\t"):
    a, b = _acc(r["Query_file"]), _acc(r["Ref_file"])
    other = a if "MF6" in (r["Ref_file"] or "") + (r["Ref_name"] or "") else b
    if "MF6" in str(other):
        continue
    ani[other] = max(ani.get(other, 0.0), float(r["ANI"]))

org, host = {}, {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    org[r["accession"]] = r["organism_name"]
    host[r["accession"]] = r["host_category"]

c3 = {r["accession"]: (r["c3_present"] == "True")
      for r in csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"),
                              delimiter="\t")}
c3["MF6"] = True

carriers = defaultdict(set)
for r in csv.DictReader(open(TAB / "toxin12_carriers.tsv"), delimiter="\t"):
    if r["tier1_carrier"] == "True":
        carriers[r["locus"]].add(r["accession"])

loci = {r["locus"]: r for r in
        csv.DictReader(open(TAB / "toxin12_locus_summary.tsv"), delimiter="\t")}

# same inner->outer ordering rule as fig11: fewest carriers first
ORDER = sorted(RING_LOCI, key=lambda k: len(carriers[k]))
prot_col = {k: PROT_COLOR[i % len(PROT_COLOR)] for i, k in enumerate(ORDER)}

# ------------------------------------------------------------------ prune
tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
# Same rooting as figures 4, 10 and 11, from the shared conventions module, so
# every figure drawn off this tree is oriented the same way. It makes no
# difference to what is PRUNED here -- the selected clade sits deep inside the
# ingroup and both candidate roots lie outside it, so the MRCA is identical
# either way -- but a figure that quietly used a different rooting from its
# neighbours is a trap for whoever reads the scripts next.
tree.root_at_midpoint()          # Bio.Phylo needs a root to walk from
_c3 = {}
for _r in csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t"):
    _c3[_r["accession"]] = (_r["c3_present"] == "True")
_c3["MF6"] = True
td.ingroup_root(tree, _c3, org, verbose=False)
tree.ladderize()
alltips = {t.name: t for t in tree.get_terminals()}

sel = [a for a, v in ani.items() if v >= ANI_MIN and a in alltips] + ["MF6"]
clade = tree.common_ancestor([alltips[x] for x in sel])
tips = clade.get_terminals()
names = [t.name for t in tips]
extra = [n for n in names if n not in sel]
print(f"selected at >={ANI_MIN}% ANI and on the tree: {len(sel)}")
print(f"MRCA clade: {len(names)} tips ({len(extra)} intervening: {extra or 'none'})")

missing = sorted((a for a, v in ani.items() if v >= ANI_MIN and a not in alltips),
                 key=lambda a: -ani[a])
print(f"relatives >={ANI_MIN}% NOT on the tree: {len(missing)} (the rebuild does not dereplicate; absentees were dropped from the core alignment)")
for a in missing[:4]:
    print(f"   {a}  ANI {ani[a]:.2f}  {org.get(a, '?')}")

sup = [c.confidence for c in clade.get_nonterminals() if c.confidence is not None]
print(f"internal nodes: {len(clade.get_nonterminals())}, support range "
      f"{min(sup):.0f}-{max(sup):.0f}")

# ------------------------------------------------------------------ layout
# rectangular phylogram: x = root-to-node branch-length distance within the
# pruned clade, y = tip order. Distances are the ORIGINAL ones -- pruning here
# only drops sisters outside the clade, it never re-estimates anything.
xpos, ypos = {}, {}
xpos[id(clade)] = 0.0
for cl in clade.find_clades(order="preorder"):
    for c in cl.clades:
        xpos[id(c)] = xpos[id(cl)] + max(c.branch_length or 0.0, 0.0)

# ---- MF7, GRAFTED -- not inferred -----------------------------------------
# MF7 is absent from the chromosome-1 alignment (D13: its 8-contig chr1 was
# typed to a single 216 kb fragment, giving 93.96% gaps), so it cannot be pruned
# out of this tree like every other tip here. It is nonetheless the second
# sequencing of MF6's own strain at ANI 100.00 (D10) and belongs in a panel
# about MF6's closest relatives -- and unlike a blank placeholder it has real
# data in every column of the matrix, which is D18's point: with its complete
# pC3 it is indistinguishable from MF6.
#
# It is therefore drawn at ZERO divergence from MF6 with a DASHED connector and
# an open star, and every statistic below is still computed over `tips` -- the
# alignment-based set -- never over the grafted row. Placement asserted from ANI,
# not from the alignment; the caption says so.
mf7 = Clade(name="MF7", branch_length=0.0)
_mf6 = next((t for t in tips if t.name == "MF6"), None)
if _mf6 is None:
    sys.exit("MF6 is not in the pruned clade - cannot graft MF7 beside it")
draw_tips = []
for t in tips:
    draw_tips.append(t)
    if t is _mf6:
        draw_tips.append(mf7)
xpos[id(mf7)] = xpos[id(_mf6)]
GRAFTED = {id(mf7)}
print(f"grafted MF7 beside MF6 at zero divergence (ANI {ani.get('MF7', float('nan')):.2f})")

for i, t in enumerate(draw_tips):
    ypos[id(t)] = float(i)
for cl in clade.find_clades(order="postorder"):
    if cl.clades:
        ypos[id(cl)] = float(np.mean([ypos[id(c)] for c in cl.clades]))
XMAX = max(xpos[id(t)] for t in draw_tips)

fig = plt.figure(figsize=(8.6, 6.4))
gsp = fig.add_gridspec(1, 2, width_ratios=[1.30, 1.0],
                       left=0.035, right=0.985, top=0.868, bottom=0.150,
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

# Support: every internal node here is 98-100, so printing 17 near-identical
# numbers would be pure clutter. Only a node below 100 gets a label; the rest
# are covered by the note in the subtitle.
for cl in clade.get_nonterminals():
    if cl.confidence is not None and cl.confidence < 100 and cl is not clade:
        axT.text(xpos[id(cl)] - XMAX * 0.006, ypos[id(cl)] + 0.30,
                 f"{cl.confidence:.0f}", fontsize=5.4, ha="right", va="center",
                 color="#777777", zorder=5)

LX = XMAX * 1.035
# dashed connector from MF6's tip down to the grafted MF7 row
axT.plot([xpos[id(mf7)]] * 2, [ypos[id(_mf6)], ypos[id(mf7)]], lw=0.9,
         ls=(0, (2, 2)), color=STAR, zorder=3)
for t in draw_tips:
    y = ypos[id(t)]
    grafted = id(t) in GRAFTED
    axT.plot([xpos[id(t)], LX - XMAX * 0.008], [y, y], lw=0.4,
             ls=(0, (2, 2)) if grafted else (0, (1, 2)),
             color=STAR if grafted else "#C8C8C8", zorder=1)
    is_mf6 = t.name == "MF6"
    sp = org.get(t.name, t.name)
    sp = (sp.replace("Burkholderia ", "B. ").replace("Paraburkholderia ", "P. "))
    if grafted:
        lab = "B. sola MF7   grafted, not in the alignment"
    elif is_mf6:
        lab = "B. sola MF6"
    else:
        lab = f"{sp}   {t.name}"
    axT.text(LX, y, lab, fontsize=6.3, va="center", ha="left",
             color=STAR if (is_mf6 or grafted) else "#222222",
             fontweight="bold" if (is_mf6 or grafted) else "normal",
             fontstyle="italic" if grafted else "normal", zorder=5)
    if is_mf6 or grafted:
        # filled star = a tip the alignment placed; OPEN star = grafted by ANI
        axT.scatter([LX - XMAX * 0.030], [y], s=52, marker="*",
                    c="none" if grafted else STAR, edgecolors=STAR if grafted else "black",
                    linewidths=0.7 if grafted else 0.4, zorder=6)
    a = ani.get(t.name)
    if a is not None:
        axT.text(LX + XMAX * 0.86, y, f"{a:.2f}", fontsize=5.8, va="center",
                 ha="right", color="#666666", zorder=5)

axT.text(LX + XMAX * 0.86, -1.15, "ANI to\nMF6 (%)", fontsize=5.8, ha="right",
         va="center", color="#666666", linespacing=1.3)
axT.set_xlim(-XMAX * 0.02, XMAX * 1.92)
# headroom above row 0 is for the rotated column headers, which are drawn in
# axM's data coordinates and would otherwise run off the top into the subtitle
axT.set_ylim(len(draw_tips) - 0.4, -5.1)
axT.axis("off")

# scale bar
_n = 10 ** np.floor(np.log10(XMAX / 3))
step = _n * max(1, round((XMAX / 3) / _n))
axT.plot([0, step], [len(draw_tips) - 0.75] * 2, lw=1.2, color="#000000")
axT.text(step / 2, len(draw_tips) - 1.05, f"{step:g} subs/site", ha="center",
         va="bottom", fontsize=6)

# ------------------------------------------------------------------ matrix
COLS = ([("pC3", None, figstyle.C3_COLORS[True]),
         ("isolation source", None, "#BBBBBB")]
        + [(ring_label(k), k, prot_col[k]) for k in ORDER])
for j, (lab, key, col) in enumerate(COLS):
    for t in draw_tips:
        y = ypos[id(t)]
        hatch = None
        if key is None and lab == "pC3":
            fc = figstyle.C3_COLORS[True] if c3.get(t.name) else ABSENT
        elif key is None:
            h = host.get(t.name, "unknown")
            fc = figstyle.HOST_COLORS.get(h, ABSENT)
            # figstyle paints "unknown" #EEEEEE, one shade off the ABSENT grey.
            # Side by side in one matrix that is unreadable -- a genome with no
            # recorded source would look like a missing value. Hatch it instead;
            # every genome HAS a source, we just do not always know it.
            if h == "unknown":
                hatch = "///"
        else:
            fc = col if t.name in carriers[key] else ABSENT
        _g = id(t) in GRAFTED
        axM.add_patch(Rectangle((j + 0.06, y - 0.42), 0.88, 0.84,
                                facecolor=fc, hatch=hatch,
                                edgecolor=STAR if _g else ("#B0B0B0" if hatch else "none"),
                                linewidth=0.6 if _g else 0.3))
    axM.text(j + 0.5, -0.75, lab, fontsize=5.9, rotation=90, ha="center",
             va="bottom", color="#222222")
# separator between the genome-property columns and the toxin block
axM.plot([2.0, 2.0], [-0.55, len(draw_tips) - 0.45], color="#FFFFFF", lw=2.4,
         zorder=4)
axM.set_xlim(0, len(COLS))
axM.axis("off")

# ------------------------------------------------------------------ labels
n_pc3 = sum(1 for t in tips if c3.get(t.name))
n_loci = {t.name: sum(1 for k in RING_LOCI if t.name in carriers[k]) for t in tips}
full = sorted(n for n in n_loci if n_loci[n] == len(RING_LOCI))
print(f"carrying all {len(RING_LOCI)} loci: {full}")
print("loci per genome: " + ", ".join(f"{n_loci[t.name]}" for t in tips))

# The first draft of this title said "every member carries at least one of the
# ten loci". True, but empty -- three of the ten are in >80% of all 772 genomes,
# so that statement is guaranteed before you look. The real result is the
# gradient: the complete set stops dead at MF6's nearest neighbour.
fig.text(0.035, 0.982,
         "MF6's closest sequenced relative carries the same complete ten-locus "
         "set;\nthe rest of the clade carries only the near-universal subset",
         fontsize=10.5, fontweight="bold", va="top", linespacing=1.35)
fig.text(0.035, 0.905,
         f"Pruned from the {len(alltips)}-tip chromosome-1 tree — branch lengths and support are as inferred there, not re-estimated. Tips are MF6 plus\n"
         f"every genome on that tree at ≥{ANI_MIN:g}% ANI to MF6; "
         + ("that set is exactly monophyletic, so no intervening lineage is hidden. "
            if not extra else
            f"{len(extra)} intervening lineage(s) fall inside the clade and are shown. ")
         + f"Internal-node\nsupport spans {min(sup):.0f}–{max(sup):.0f}% UFBoot across "
         f"{len(clade.get_nonterminals())} nodes (only values below 100 are printed). "
         f"pC3 is present in {n_pc3}/{len(tips)}. Loci per genome:\n"
         + ", ".join(str(n) for n in sorted((n_loci[t.name] for t in tips), reverse=True)) + ".",
         fontsize=6.3, color="#444444", va="top", linespacing=1.5)

# The original hard-coded exactly two absentees, both lost to dereplication.
# The rebuild dereplicates nothing, so whoever is missing is missing for a
# different reason and the count is no longer fixed -- build the note from the
# data or omit it entirely.
# MF7 is now GRAFTED into the panel rather than listed as absent, so it must not
# also be announced as "not shown" -- that note contradicted the figure the
# moment the graft was added. Anyone else still missing keeps the old wording.
_still_missing = [a for a in missing if a != "MF7"]
if "MF7" in missing:
    fig.text(0.035, 0.088,
             f"MF7 is GRAFTED, not inferred (open star, dashed): it is absent from "
             f"the chromosome-1 alignment — D13 dropped it at 93.96% gaps, its "
             f"8-contig\nchromosome 1 having been typed to a single 216 kb fragment "
             f"— so there is nothing to prune to. It is drawn at zero divergence "
             f"from MF6 on the\nstrength of ANI 100.00 (D10), which is a placement "
             f"asserted from ANI and not from this alignment. Its matrix row is "
             f"measured data, and it\nmatches MF6 across all ten loci (D18). Every "
             f"statistic quoted above is computed over the {len(tips)} inferred "
             f"tips only, never the grafted row.",
             fontsize=5.9, color=STAR, va="top", linespacing=1.5)
elif _still_missing:
    _who = ", ".join(f"{a} ({ani[a]:.2f}% ANI)" for a in _still_missing[:3])
    fig.text(0.035, 0.088,
             f"NOT SHOWN: {_who} — ≥{ANI_MIN:g}% ANI to MF6 but absent from the "
             f"chromosome-1 tree, having been dropped from the core alignment\n"
             f"(D13: >90% gaps), so there is nothing to prune to. Still present in "
             f"Figure 11 panel B. Placing them here needs a fresh inference from "
             f"genomes, which this figure does not do.",
             fontsize=5.9, color=STAR, va="top", linespacing=1.5)

hosts_here = [h for h in figstyle.HOST_COLORS
              if any(host.get(t.name) == h for t in tips)]
leg = fig.legend(handles=[Patch(facecolor=figstyle.C3_COLORS[True],
                                label="pC3 present"),
                          Patch(facecolor=ABSENT, label="absent (any column)")]
                 + [Patch(facecolor=figstyle.HOST_COLORS[h],
                          edgecolor="#B0B0B0" if h == "unknown" else "#999999",
                          hatch="///" if h == "unknown" else None, linewidth=0.3,
                          label="source not recorded" if h == "unknown"
                          else h.replace("_", " ")) for h in hosts_here],
                 loc="lower center", ncol=8, frameon=False, fontsize=6.2,
                 handlelength=1.2, columnspacing=1.1,
                 bbox_to_anchor=(0.5, -0.004),
                 title="isolation source · toxin columns use the Figure 11 "
                       "locus colours",
                 title_fontsize=6.8)
leg.get_title().set_fontweight("bold")

ok = figstyle.save(fig, str(FIG / "fig12_close_relatives"))
print("OK" if ok else "FAILED vector-text check")
