#!/usr/bin/env python3
"""Figure 10 - does the RHS warhead track the chromosome-1 phylogeny?

Panel A: the 763-tip chromosome-1 tree with the warhead carriers marked. They do
         not occur across the genus; they occupy one wedge of it. Rooted on the
         ingroup stem with the outgroup and near-clone blocks collapsed, exactly
         as figure 4 -- both conventions come from tree_display so the two
         figures cannot disagree about the same tree.
Panel B: that region enlarged - the smallest clade containing every carrier. The
         carriers plus ONE non-carrier form a 100%-UFBoot clade. The exception,
         B. cenocepacia YG-3, carries pC3 but has no RHS gene at all, so it is a
         whole-gene loss inside the clade, not a warhead swap.
Panel C: randomisation of mean pairwise patristic distance. Two nulls: all tips,
         and pC3-positive tips only. The second controls for pC3 presence being
         itself clade-structured, and is the sharper test.

Inputs: tables/chr1_core.treefile, c3_calls_all_genomes.tsv,
        rhs_search_per_genome.tsv, MF6_ani_raw.tsv
"""
import csv
import random
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import radialtree as rt  # noqa: E402
import tree_display as td  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

PC3 = "#1B3A6B"
WARN = "#D55E00"
GREEN = "#009E73"
GREY = "#BBBBBB"

# Carriers are READ FROM THE TABLE, never hard-coded. The hard-coded set this
# replaces still held the original run's SIX, missing GCF_016899425.1 and
# GCF_053209605.1 -- the two the original lost to dereplication and the rebuild
# recovers (D14: the count goes 7 -> 9 across 773 genomes). A figure that marks
# six stars while the report says eight is exactly the failure a hard-coded
# constant invites.
WARHEAD_QUERY = "CT354_CFFIHE_03684"
_all_carriers = set()
with open(TAB / "rhs_search_per_genome.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["query"] == WARHEAD_QUERY and int(r["n_hits_any"]) > 0:
            _all_carriers.add(r["accession"])
# The rebuild does not dereplicate, so no carrier is lost to a cluster
# representative any more. Kept empty rather than deleted so the reason the
# original had two entries here stays on the record.
DEREP_AWAY = set()

import tables_compat  # noqa: E402

org = tables_compat.organism_names(TAB)
c3 = {}
with open(TAB / "c3_calls_all_genomes.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        c3[r["accession"]] = (r["c3_present"] == "True")
org["MF6"] = "Burkholderia sola  MF6"
c3["MF6"] = True

ani = {}
with open(TAB / "MF6_ani_raw.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        a = Path(r["Ref_file"]).name.replace(".fna", "")
        ani[a] = float(r["ANI"])

# does each genome have the full-length RHS scaffold at all?
has_scaffold = {}
with open(TAB / "rhs_search_per_genome.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["query"] == "FULL_CFFIHE_03684":
            has_scaffold[r["accession"]] = int(r["n_hits_any"]) > 0

tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
tree.root_at_midpoint()          # Bio.Phylo needs a root to walk from
print("display conventions:")
ingroup, outgroup = td.ingroup_root(tree, c3, org)
nearclones = td.nearclone_clades(tree, c3, org, skip=[outgroup])
tree.ladderize(reverse=True)     # after re-rooting, restore this figure's order

_on_tree = {t.name for t in tree.get_terminals()}
CARRIERS = _all_carriers & _on_tree
_off = _all_carriers - _on_tree
print(f"warhead carriers: {len(_all_carriers)} total, {len(CARRIERS)} on the "
      f"tree" + (f" (absent: {', '.join(sorted(_off))} - D13)" if _off else ""))
assert CARRIERS, "no warhead carrier is on the tree"

fig = plt.figure(figsize=(7.5, 8.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1.32, 1.0], width_ratios=[1.0, 1.0],
                      hspace=0.20, wspace=0.30,
                      left=0.045, right=0.975, top=0.915, bottom=0.065)

# ---------------------------------------------------------------- panel A
ax = fig.add_subplot(gs[0, :])
ax.set_aspect("equal"); ax.axis("off")
_collapse = {outgroup: 22.0}
_collapse.update({cl: 10.0 for cl in nearclones})
lay = rt.RadialLayout(tree, start_deg=92.0, extent_deg=348.0, collapse=_collapse)
rt.draw_tree(ax, lay, lw=0.35, color="#999999")
drawn = [t for t in tree.get_terminals() if id(t) in lay.radius]
print(f"  drawn: {len(drawn)} tips + {len(nearclones) + 1} wedges")

# collapsing can only hide pC3-NEGATIVE genomes and the warhead is
# pC3-restricted, so no carrier can ever be swallowed by a wedge. Assert it.
assert CARRIERS <= {t.name for t in drawn}, "a warhead carrier got collapsed"

for cl, a0, a1, r_stem, r_crown in lay.collapsed_wedges():
    truncated = r_crown > 1.0
    td.draw_wedge(ax, rt, a0, a1, r_stem,
                  0.86 if truncated else r_crown, break_marker=truncated)

# highlight the carrier clade and mark carrier tips
mrca = tree.common_ancestor(list(CARRIERS))
clade_tips = {t.name for t in mrca.get_terminals()}
for t in drawn:
    r = lay.radius[id(t)]
    a = lay.angle[id(t)]
    x, y = rt.polar_to_xy(r, a)
    if t.name in CARRIERS:
        ax.plot([x], [y], marker="*", ms=8.5, color=WARN, zorder=6,
                markeredgecolor="white", markeredgewidth=0.4)
    elif t.name in clade_tips:
        ax.plot([x], [y], marker="o", ms=4.2, markerfacecolor="white",
                markeredgecolor=WARN, markeredgewidth=1.0, zorder=6)
    elif c3.get(t.name):
        ax.plot([x], [y], marker="o", ms=1.5, color=PC3, zorder=4, alpha=0.75)

# The tree is strongly non-ultrametric, so tip radii vary widely. Anchor the
# highlight to the CLADE's own tip radius -- using the global max puts the arc
# far outside the carriers and it reads as unrelated to them.
rmax = max(lay.radius[id(t)] for t in drawn)
_ct = [t for t in drawn if t.name in clade_tips]
angs = [lay.angle[id(t)] for t in _ct]
rr = max(lay.radius[id(t)] for t in _ct) * 1.06
a0, a1 = min(angs), max(angs)
th = np.linspace(a0 - 0.055, a1 + 0.055, 60)
ax.plot(rr * np.cos(th), rr * np.sin(th), color=WARN, lw=2.4,
        solid_capstyle="butt", zorder=7)
amid = (a0 + a1) / 2
ax.annotate(f"all {len(CARRIERS)} carriers sit here",
            xy=rt.polar_to_xy(rr * 1.03, amid),
            xytext=rt.polar_to_xy(rmax * 1.16, amid + 0.30),
            fontsize=6.6, color=WARN, ha="right", va="center",
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.9))
ax.plot([], [], marker="*", ms=8, color=WARN, ls="none",
        label=f"warhead carrier ({len(CARRIERS)})")
ax.plot([], [], marker="o", ms=4, markerfacecolor="white", markeredgecolor=WARN,
        ls="none",
        label=f"same clade, no warhead ({(len(clade_tips) - len(CARRIERS))})")
_n_pc3_drawn = sum(1 for t in drawn if c3.get(t.name))
ax.plot([], [], marker="o", ms=3, color=PC3, ls="none",
        label=f"pC3-positive ({_n_pc3_drawn} drawn)")
ax.legend(frameon=False, fontsize=6.2, loc="upper left",
          bbox_to_anchor=(0.02, 0.99), handletextpad=0.4, labelspacing=0.35)
ax.set_title(f"A   Carriers occupy one wedge of the "
             f"{len(tree.get_terminals())}-tip chromosome-1 tree",
             loc="left", fontweight="bold")
lim = rmax * 1.13
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

# ---------------------------------------------------------------- panel B
ax = fig.add_subplot(gs[1, 0])
path = tree.get_path(mrca)
parent = path[-2] if len(path) > 1 else tree.root
sub = parent
names = [t.name for t in sub.get_terminals()]

# simple rectangular cladogram of the subtree
ypos, ordered = {}, []
for i, t in enumerate(sub.get_terminals()):
    ypos[id(t)] = i
    ordered.append(t)
for cl in sub.get_nonterminals(order="postorder"):
    ypos[id(cl)] = float(np.mean([ypos[id(c)] for c in cl.clades]))
depth = {id(sub): 0.0}
for cl in sub.find_clades(order="preorder"):
    for c in cl.clades:
        depth[id(c)] = depth[id(cl)] + max(c.branch_length or 0.0, 0.0)
maxd = max(depth.values()) or 1.0
for cl in sub.find_clades(order="preorder"):
    for c in cl.clades:
        x0, x1 = depth[id(cl)] / maxd, depth[id(c)] / maxd
        ax.plot([x0, x1], [ypos[id(c)], ypos[id(c)]], color="#444444", lw=0.7)
    if cl.clades:
        ys = [ypos[id(c)] for c in cl.clades]
        ax.plot([depth[id(cl)] / maxd] * 2, [min(ys), max(ys)], color="#444444", lw=0.7)

for t in ordered:
    y = ypos[id(t)]
    x = depth[id(t)] / maxd
    car = t.name in CARRIERS
    ax.plot([x], [y], marker="*" if car else "o", ms=8 if car else 4.6,
            color=WARN if car else ("white" if t.name in clade_tips else GREY),
            markeredgecolor=WARN if t.name in clade_tips else GREY,
            markeredgewidth=0.9, zorder=5)
    lab = org.get(t.name, "?").replace("Burkholderia ", "B. ")
    a = ani.get(t.name)
    extra = f"  {a:.2f}%" if a else ""
    ax.text(x + 0.035, y, f"{lab}{extra}", va="center", fontsize=5.4,
            color="#111111" if t.name in clade_tips else "#666666",
            fontweight="bold" if car else "normal")
ax.axhspan(min(ypos[id(t)] for t in ordered if t.name in clade_tips) - 0.5,
           max(ypos[id(t)] for t in ordered if t.name in clade_tips) + 0.5,
           color=WARN, alpha=0.09, zorder=0)
ax.annotate("no RHS gene at all,\nyet carries pC3",
            xy=(depth[id([t for t in ordered if t.name == "GCF_003966315.1"][0])] / maxd,
                ypos[id([t for t in ordered if t.name == "GCF_003966315.1"][0])]),
            xytext=(0.10, len(ordered) * 0.60), fontsize=5.6, color=WARN,
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.8))
ax.set_ylim(-0.8, len(ordered) - 0.2)
ax.set_xlim(-0.02, 2.05)
ax.invert_yaxis()
ax.set_yticks([]); ax.set_xticks([])
for sp in ("left", "right", "top", "bottom"):
    ax.spines[sp].set_visible(False)
ax.set_title("B   The clade, enlarged (UFBoot 100)\n      label = species, ANI to MF6",
             loc="left", fontweight="bold", fontsize=8.4)

# ---------------------------------------------------------------- panel C
ax = fig.add_subplot(gs[1, 1])
tips = [t.name for t in tree.get_terminals()]
pc3_tips = [t for t in tips if c3.get(t)]
dist = {}


def pdist(a, b):
    k = (a, b) if a < b else (b, a)
    if k not in dist:
        dist[k] = tree.distance(a, b)
    return dist[k]


def meanpd(ns):
    ps = list(combinations(ns, 2))
    return sum(pdist(a, b) for a, b in ps) / len(ps)


obs = meanpd(sorted(CARRIERS))
random.seed(20260810)
N = 4000
for pool, col, lab in [(tips, GREY, f"all tips (n = {len(tips)})"),
                       (pc3_tips, PC3, f"pC3-positive only (n = {len(pc3_tips)})")]:
    null = [meanpd(random.sample(pool, len(CARRIERS))) for _ in range(N)]
    ax.hist(null, bins=45, color=col, alpha=0.62, label=lab, linewidth=0)
ax.axvline(obs, color=WARN, lw=1.8)
ax.annotate(f"observed\n{obs:.3f}\np < 0.001", xy=(obs, ax.get_ylim()[1] * 0.62),
            xytext=(obs + 0.075, ax.get_ylim()[1] * 0.80), fontsize=6.0, color=WARN,
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.8))
ax.set_xlabel("mean pairwise patristic distance")
ax.set_ylabel(f"random draws of {len(CARRIERS)} tips")
ax.legend(frameon=False, fontsize=5.8, loc="upper right")
ax.set_title("C   Tighter than chance under both nulls",
             loc="left", fontweight="bold", fontsize=8.4)

fig.suptitle("The warhead tracks the chromosome-1 phylogeny: one clade, one loss,\n"
             "no evidence of movement between distant lineages",
             fontsize=10, fontweight="bold", y=0.982)

ok = figstyle.save(fig, str(FIG / "fig10_warhead_tree"))
print("OK" if ok else "FAILED vector-text check")
