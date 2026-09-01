#!/usr/bin/env python3
"""Figure - the c3 (pC3) phylogeny, annotated with host and replicon size.

Radial maximum-likelihood phylogram of the 140 replicons that carry a genuine
pC3 gene-content signature, built by IQ-TREE from the 45-family c3 core-gene
alignment (GTR+F+I+R5, 1000 ultrafast bootstraps).

Rings, inner to outer:
  1. species group   -- B. cepacia complex vs other Burkholderia vs other genera
  2. host category   -- from NCBI BioSample, auto-curated

The c3-length bars were removed 2026-08-09 at Moshe's request. Dotted leaders
now run from each tip out to the rings: this is a phylogram, so tips end at
different radii and without them the rings float free of the tree.

Only c3-BEARING genomes are on this tree, so presence/absence cannot be drawn
here; that mapping needs the chromosome-1 species tree, which is still running.

Node dots mark UFBoot >= 95.
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import radialtree as rt  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

TREE = TAB / "c3_core.treefile"
if not TREE.exists():
    sys.exit(f"{TREE} not found - run pull_results.sh")

# ---------------------------------------------------------------- metadata
hosts = {r["accession"]: r for r in
         csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t")}
# MF6 is a lab isolate with no NCBI BioSample; it is carried as a manually
# curated row in host_categories.tsv so this figure and fig1 agree. Before that
# row existed MF6 fell through to "unknown" here while fig1 hard-coded it.
if "MF6" not in hosts:
    print("WARNING: MF6 missing from host_categories.tsv - it will draw as unknown")
else:
    print(f"MF6 host: {hosts['MF6']['host_category']}")

# B. cepacia complex sensu lato -- the clade Agnoli et al. surveyed
BCC = {"cepacia", "multivorans", "cenocepacia", "stabilis", "vietnamiensis",
       "dolosa", "ambifaria", "anthina", "pyrrocinia", "ubonensis", "latens",
       "diffusa", "arboris", "seminalis", "metallica", "contaminans", "lata",
       "pseudomultivorans", "stagnalis", "territorii", "puraquae", "orbicola",
       "aenigmatica", "sola", "semiarida", "catudaia", "alpina"}


def group(acc):
    if acc == "MF6":
        return "bcc"          # established as B. sola by ANI
    o = hosts.get(acc, {}).get("organism_name", "")
    parts = o.split()
    if not parts:
        return "other_genus"
    gen, sp = parts[0], (parts[1] if len(parts) > 1 else "")
    if gen != "Burkholderia":
        return "other_genus"
    return "bcc" if sp in BCC else "other_burk"


GRP_COLORS = {"bcc": "#1B3A6B", "other_burk": "#7FA6D9", "other_genus": "#D9D9D9"}
GRP_LABEL = {"bcc": "B. cepacia complex", "other_burk": "other Burkholderia",
             "other_genus": "other genus"}

# ---------------------------------------------------------------- tree
tree = Phylo.read(str(TREE), "newick")
tree.root_at_midpoint()
tree.ladderize(reverse=True)
tips = tree.get_terminals()
print(f"tips: {len(tips)}")

lay = rt.RadialLayout(tree, start_deg=92.0, extent_deg=348.0)

fig = plt.figure(figsize=(7.2, 7.2))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.92])
ax.set_aspect("equal")
ax.axis("off")

rt.draw_tree(ax, lay, lw=0.4, color="#444444")

# UFBoot >= 95 as small dots. IQ-TREE writes support in the node's confidence
# field for a newick with internal labels; Bio.Phylo puts it in .confidence.
n_sup = 0
for cl in tree.get_nonterminals():
    c = cl.confidence
    if c is None:
        continue
    if c >= 95:
        x, y = rt.polar_to_xy(lay.radius[id(cl)], lay.angle[id(cl)])
        ax.plot([x], [y], marker="o", ms=1.5, mfc="#000000", mec="none", zorder=4)
        n_sup += 1
print(f"nodes with UFBoot >= 95: {n_sup}")

R0 = 1.02

# dotted leaders from each tip out to the rings
_ext = [[rt.polar_to_xy(lay.tip_radius(t), lay.tip_angle(t)),
         rt.polar_to_xy(R0, lay.tip_angle(t))] for t in tips]
ax.add_collection(LineCollection(_ext, colors="#BBBBBB", linewidths=0.22,
                                 linestyles=(0, (1, 2)), zorder=1))

rt.draw_ring(ax, lay, {t.name: group(t.name) for t in tips},
             GRP_COLORS, R0, R0 + 0.05)
rt.draw_ring(ax, lay,
             {t.name: hosts.get(t.name, {}).get("host_category", "unknown")
              for t in tips},
             figstyle.HOST_COLORS, R0 + 0.06, R0 + 0.11)

# MF6 callout
mf6 = next((t for t in tips if t.name == "MF6"), None)
if mf6 is not None:
    a = lay.tip_angle(mf6)
    x0, y0 = rt.polar_to_xy(R0 + 0.13, a)
    x1, y1 = rt.polar_to_xy(R0 + 0.22, a)
    ax.plot([x0, x1], [y0, y1], lw=0.9, color="#D55E00", zorder=5)
    ax.text(*rt.polar_to_xy(R0 + 0.30, a), "MF6\n(B. sola)", ha="center",
            va="center", fontsize=7, color="#D55E00", fontweight="bold",
            zorder=5, style="italic")

LIM = 1.42
ax.set_xlim(-LIM, LIM)
ax.set_ylim(-LIM, LIM)

# scale bar in substitutions/site (radii were normalised to the tree depth)
depth = lay.max_r
nice = 10 ** np.floor(np.log10(depth / 4))
step = nice * max(1, round((depth / 4) / nice))
ax.plot([-LIM + 0.04, -LIM + 0.04 + step / depth], [-LIM + 0.06, -LIM + 0.06],
        lw=1.2, color="#000000")
ax.text(-LIM + 0.04 + step / (2 * depth), -LIM + 0.11, f"{step:g} subs/site",
        ha="center", va="bottom", fontsize=7)

# only groups actually on the tree -- no genus outside Burkholderia sensu
# stricto carries a c3, so "other genus" would be an empty legend entry
_on_tree = {group(t.name) for t in tips}
leg1 = [Patch(facecolor=GRP_COLORS[k], label=GRP_LABEL[k])
        for k in ("bcc", "other_burk", "other_genus") if k in _on_tree]
present = [h for h in figstyle.HOST_COLORS
           if any(hosts.get(t.name, {}).get("host_category") == h for t in tips)]
leg2 = [Patch(facecolor=figstyle.HOST_COLORS[h], label=h.replace("_", " "))
        for h in present]
leg3 = [Line2D([], [], marker="o", ls="none", ms=3, color="#000000",
               label="UFBoot $\\geq$ 95")]

l1 = ax.legend(handles=leg1, loc="upper left", bbox_to_anchor=(-0.01, 1.005),
               frameon=False, fontsize=6.8, title="Species group",
               title_fontsize=7, handlelength=1.0, handleheight=0.9)
l1.get_title().set_fontweight("bold")
l2 = ax.legend(handles=leg2, loc="upper right", bbox_to_anchor=(1.01, 1.005),
               frameon=False, fontsize=6.8, title="Host / source",
               title_fontsize=7, handlelength=1.0, handleheight=0.9)
l2.get_title().set_fontweight("bold")
ax.add_artist(l1)
ax.legend(handles=leg3, loc="lower right", bbox_to_anchor=(1.02, -0.01),
          frameon=False, fontsize=6.8, handlelength=1.0, handleheight=0.9)
ax.add_artist(l2)

fig.suptitle("c3 (pC3) core-gene phylogeny of 140 Burkholderia sensu lato genomes",
             fontsize=10, fontweight="bold", y=0.985)

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig6_c3_phylogeny"))
print("done")
