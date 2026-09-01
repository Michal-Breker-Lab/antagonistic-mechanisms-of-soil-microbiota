#!/usr/bin/env python3
"""Figure 3 - species tree vs c3 tree.

Congruent topologies mean c3 has been inherited vertically with its host
chromosome; incongruence means it has moved horizontally between lineages. This
comparison is the substantive reason for building two trees rather than one.

Connector lines are coloured by host category so horizontal movement between
habitats, if present, is visible directly.

The right-hand tree is rotated to reduce crossings before drawing. Rotation
changes only the drawing, never the topology -- swapping the two children of a
node yields the same tree. Without it a 300-tip tanglegram is unreadable, and
the apparent number of crossings would be an artifact of arbitrary node order.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

T1 = TAB / "chr1_core.treefile"
T2 = TAB / "c3_core.treefile"
for t in (T1, T2):
    if not t.exists():
        sys.exit(f"missing tree: {t}")

# IQ-TREE trees are unrooted. Midpoint-root both before laying them out --
# without a root the horizontal (depth) axis is meaningless and the two trees
# cannot be visually compared. Rooting is a display choice; the Robinson-Foulds
# statistic reported in the phylogenetic-statistics table is computed on the
# unrooted topologies.
sp = Phylo.read(str(T1), "newick"); sp.root_at_midpoint(); sp.ladderize()
c3 = Phylo.read(str(T2), "newick"); c3.root_at_midpoint(); c3.ladderize()

shared = sorted({t.name for t in sp.get_terminals()} &
                {t.name for t in c3.get_terminals()})
print(f"species tree tips: {len(sp.get_terminals())}, "
      f"c3 tree tips: {len(c3.get_terminals())}, shared: {len(shared)}")
if len(shared) < 5:
    sys.exit("too few shared tips")

keep = set(shared)
for tr in (sp, c3):
    for t in [x for x in tr.get_terminals() if x.name not in keep]:
        tr.prune(t)
sp.root_at_midpoint(); sp.ladderize()
c3.root_at_midpoint(); c3.ladderize()


def tip_order(tree):
    return [t.name for t in tree.get_terminals()]


def untangle(tree, target_rank, rounds=6):
    """Rotate nodes to reduce crossings against target_rank. Topology-preserving."""
    for _ in range(rounds):
        changed = False
        for cl in tree.get_nonterminals(order="postorder"):
            if len(cl.clades) != 2:
                continue
            def mean_rank(node):
                rs = [target_rank.get(t.name) for t in node.get_terminals()]
                rs = [r for r in rs if r is not None]
                return np.mean(rs) if rs else 0.0
            a, b = cl.clades
            if mean_rank(a) > mean_rank(b):
                cl.clades = [b, a]
                changed = True
        if not changed:
            break
    return tree


rank1 = {n: i for i, n in enumerate(tip_order(sp))}
c3 = untangle(c3, rank1)
order1, order2 = tip_order(sp), tip_order(c3)
rank2 = {n: i for i, n in enumerate(order2)}

# crossings = inversions between the two orderings
seq = [rank2[n] for n in order1]
crossings = sum(1 for i in range(len(seq)) for j in range(i + 1, len(seq))
                if seq[i] > seq[j])
maxc = len(seq) * (len(seq) - 1) / 2
print(f"crossings after untangling: {crossings} / {maxc:.0f} "
      f"({100*crossings/maxc:.1f}% of maximum)")
rho = np.corrcoef(np.arange(len(order1)), seq)[0, 1]
print(f"Spearman-like rank correlation of tip orders: {rho:.3f}")

host = {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    host[r["accession"]] = r["host_category"]
host.setdefault("MF6", "rhizosphere")


def layout(tree, flip=False):
    """-> dict name->(x,y) for tips, plus segment list, x in [0,1]."""
    tips = tree.get_terminals()
    ypos = {id(t): i for i, t in enumerate(tips)}
    for cl in tree.get_nonterminals(order="postorder"):
        ypos[id(cl)] = float(np.mean([ypos[id(c)] for c in cl.clades]))
    depth = {id(tree.root): 0.0}
    for cl in tree.get_nonterminals(order="preorder"):
        for c in cl.clades:
            depth[id(c)] = depth[id(cl)] + max(c.branch_length or 0.0, 0.0)
    dmax = max(depth.values()) or 1.0
    segs = []
    for cl in tree.get_nonterminals(order="preorder"):
        x0 = depth[id(cl)] / dmax
        ys = [ypos[id(c)] for c in cl.clades]
        segs.append([(x0, min(ys)), (x0, max(ys))])
        for c in cl.clades:
            segs.append([(x0, ypos[id(c)]), (depth[id(c)] / dmax, ypos[id(c)])])
    if flip:
        segs = [[(1 - x, y) for x, y in s] for s in segs]
    tipxy = {t.name: ((1 - depth[id(t)] / dmax) if flip else depth[id(t)] / dmax,
                      ypos[id(t)]) for t in tips}
    return tipxy, segs


tip1, seg1 = layout(sp)
tip2, seg2 = layout(c3, flip=True)

n = len(shared)
fig, ax = plt.subplots(figsize=(7.2, max(5.0, min(10.5, n * 0.068))))
ax.axis("off")

GAP, WID = 0.62, 1.0
s1 = [[(x * WID, y) for x, y in s] for s in seg1]
s2 = [[(WID + GAP + x * WID, y) for x, y in s] for s in seg2]
ax.add_collection(LineCollection(s1, colors="#333333", linewidths=0.35))
ax.add_collection(LineCollection(s2, colors="#333333", linewidths=0.35))

conn, cols = [], []
for name in shared:
    x1, y1 = tip1[name]
    x2, y2 = tip2[name]
    conn.append([(x1 * WID, y1), (WID + GAP + x2 * WID, y2)])
    cols.append(figstyle.HOST_COLORS.get(host.get(name, "unknown"), "#DDDDDD"))
ax.add_collection(LineCollection(conn, colors=cols, linewidths=0.55, alpha=0.75))

for name, mark in ((n, n) for n in ["MF6"]):
    if name in tip1:
        x1, y1 = tip1[name]
        x2, y2 = tip2[name]
        ax.scatter([x1 * WID, WID + GAP + x2 * WID], [y1, y2], s=30, marker="*",
                   c="#D55E00", edgecolors="black", linewidths=0.4, zorder=6)
        ax.text(x1 * WID - 0.03, y1, "MF6", fontsize=7, fontweight="bold",
                color="#D55E00", ha="right", va="center")

ax.set_xlim(-0.30, 2 * WID + GAP + 0.10)
ax.set_ylim(-4, n + 2)
ax.text(WID / 2, n + 0.5, "Chromosome-1 species tree", ha="center",
        fontsize=9, fontweight="bold")
ax.text(WID + GAP + WID / 2, n + 0.5, "c3 tree", ha="center",
        fontsize=9, fontweight="bold")
ax.text(WID + GAP / 2, -3.0,
        f"{crossings:,} crossings ({100*crossings/maxc:.1f}% of maximum); "
        f"tip-order correlation r = {rho:.2f}",
        ha="center", fontsize=7.5, color="#444444")

present = [h for h in figstyle.HOST_COLORS if h in set(host.get(k) for k in shared)]
# figure-level legend under the axes; an axes-level one collides with the
# crossings annotation at the bottom of the plot
fig.subplots_adjust(bottom=0.085, top=0.945)
fig.legend(handles=[Patch(facecolor=figstyle.HOST_COLORS[h],
                          label=h.replace("_", " ")) for h in present],
           frameon=False, fontsize=6.4, ncol=5, loc="lower center",
           bbox_to_anchor=(0.5, 0.0), handlelength=1.2, columnspacing=1.2)
fig.suptitle("Vertical inheritance or horizontal movement? "
             "Species tree versus c3 tree",
             fontsize=10.5, fontweight="bold", y=0.985)

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig3_tanglegram"))
print("done")
