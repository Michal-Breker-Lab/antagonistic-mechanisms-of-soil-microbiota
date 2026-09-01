#!/usr/bin/env python3
"""Figure 5 - replicon architecture across Burkholderia sensu lato, MF6 marked.

Panel A: size of every replicon >=100 kb, by rank within its genome, showing the
         three-replicon architecture and where c3 sits.
Panel B: architecture class by genus.
Panel C: threshold sensitivity for the "3 large replicons" call.
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402  (must import before pyplot)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB = D / "tables"
FIG = D / "figures"
FIG.mkdir(exist_ok=True)

MF6 = [3_592_343, 2_983_592, 1_179_882, 33_003, 9_111]

import rebuild_rules

rows = list(csv.DictReader(open(TAB / "replicon_census.tsv"), delimiter="\t"))
# D10: draft assemblies have contig sizes, not replicon sizes, so they cannot
# enter a size-by-rank plot or an architecture count. See rebuild_rules.
rows, _dropped = rebuild_rules.drop_drafts(rows)
print(f"genomes: {len(rows)}  (excluded {_dropped} draft assembly/assemblies: "
      f"{sorted(rebuild_rules.DRAFT_GENOMES)})")

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.30)

# ---------------- Panel A: replicon size by rank ----------------
ax = fig.add_subplot(gs[0, :])
by_rank = collections.defaultdict(list)
for r in rows:
    for i, s in enumerate(r["replicon_sizes"].split(";"), 1):
        s = int(s)
        if s >= 100_000 and i <= 5:
            by_rank[i].append(s / 1e6)

ranks = sorted(by_rank)
parts = ax.violinplot([by_rank[k] for k in ranks], positions=ranks,
                      widths=0.75, showextrema=False, showmedians=True)
for b in parts["bodies"]:
    b.set_facecolor("#9DB4D0"); b.set_edgecolor("#3A5A85"); b.set_alpha(0.85)
    b.set_linewidth(0.6)
parts["cmedians"].set_color("#1B3A6B"); parts["cmedians"].set_linewidth(1.2)

for k in ranks:
    ax.scatter(np.random.normal(k, 0.055, len(by_rank[k])),
               by_rank[k], s=1.4, c="#33445C", alpha=0.28, linewidths=0, zorder=3)

for i, s in enumerate(MF6[:3], 1):
    ax.scatter([i], [s / 1e6], s=52, marker="*", c="#D55E00",
               edgecolors="black", linewidths=0.5, zorder=6,
               label=r"MF6 ($\it{B.\ sola}$)" if i == 1 else None)

ax.set_xticks(ranks)
ax.set_xticklabels([f"{k}" for k in ranks])
ax.set_xlabel("Replicon rank within genome (1 = largest)")
ax.set_ylabel("Replicon size (Mb)")
ax.set_title(r"A   Replicon size structure across $\bf{\it{Burkholderia}}$ "
             + r"$\bf{sensu\ lato}$ " + f"(n = {len(rows)} complete genomes)",
             loc="left", fontweight="bold")
ax.axhline(0.3, color="#B00020", lw=0.8, ls="--", zorder=1)
ax.text(5.42, 0.33, "300 kb\nthreshold", fontsize=6.5, color="#B00020", va="bottom")
ax.legend(frameon=False, loc="upper right")
for k in ranks:
    ax.text(k, max(by_rank[k]) + 0.16, f"n={len(by_rank[k])}",
            ha="center", fontsize=6.5, color="#444444")
ax.set_ylim(0, 8.2)

# ---------------- Panel B: architecture by genus ----------------
ax = fig.add_subplot(gs[1, 0])
g = collections.defaultdict(collections.Counter)
for r in rows:
    g[r["genus"]][r["architecture"]] += 1
order = sorted(g, key=lambda x: -sum(g[x].values()))
classes = ["1_large", "2_large", "3_large", "4+_large"]
cols = ["#E8E8E8", "#9DB4D0", "#1B3A6B", "#5C3A6B"]
bottom = np.zeros(len(order))
for cl, c in zip(classes, cols):
    v = np.array([100 * g[x][cl] / sum(g[x].values()) for x in order])
    ax.barh(range(len(order)), v, left=bottom, color=c, height=0.72,
            edgecolor="white", linewidth=0.5, label=cl.replace("_", " "))
    bottom += v
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f"$\\it{{{x}}}$\n(n={sum(g[x].values())})" for x in order],
                   fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel("Genomes (%)")
ax.set_xlim(0, 100)
ax.set_title("B   Architecture by genus", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, ncol=2, loc="lower center",
          bbox_to_anchor=(0.5, -0.42))

# ---------------- Panel C: threshold sensitivity ----------------
ax = fig.add_subplot(gs[1, 1])
ths = [("200kb", 200), ("300kb", 300), ("500kb", 500), ("1000kb", 1000)]
counts = {}
for key, kb in ths:
    c = collections.Counter()
    for r in rows:
        n = int(r[f"n_large_{key}"])
        c[f"{n}_large" if n <= 3 else "4+_large"] += 1
    counts[kb] = c
x = np.arange(len(ths))
bottom = np.zeros(len(ths))
for cl, c in zip(classes, cols):
    v = np.array([counts[kb][cl] for _, kb in ths])
    ax.bar(x, v, bottom=bottom, color=c, width=0.68,
           edgecolor="white", linewidth=0.5)
    bottom += v
ax.set_xticks(x); ax.set_xticklabels([t[0] for t in ths])
ax.set_xlabel("Size threshold for a 'large' replicon")
ax.set_ylabel("Genomes")
ax.set_title("C   Threshold sensitivity", loc="left", fontweight="bold")
for i, (_, kb) in enumerate(ths):
    ax.text(i, counts[kb]["3_large"] / 2 + counts[kb]["1_large"] + counts[kb]["2_large"],
            str(counts[kb]["3_large"]), ha="center", va="center",
            fontsize=7, color="white", fontweight="bold")

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig5_replicon_architecture"))
print("done")
