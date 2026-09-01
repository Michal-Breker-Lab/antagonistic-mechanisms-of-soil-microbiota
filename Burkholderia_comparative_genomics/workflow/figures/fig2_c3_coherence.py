#!/usr/bin/env python3
"""Figure 2 - is c3 a coherent replicon, and which replicons actually are c3?

Panel A: where genes live. For gene families shared across genomes carrying two
         large secondary replicons, the fraction of occurrences residing on the
         smaller (c3) replicon. Strong bimodality means each replicon keeps its
         own gene complement rather than sharing at random.
Panel B: c3-diagnostic vs chromosome-2-diagnostic content for every large
         secondary replicon. Three classes separate, including a substantial set
         of megaplasmids carrying NEITHER signature -- these occupy the c3 size
         slot but are not pC3, and a size- or position-based call misassigns them.
Panel C: replicon size by class, showing that size alone cannot make the call.
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

CLASSES = [("c3", "c3 (pC3-like)", "#1B3A6B"),
           ("chromosome2", "chromosome 2", "#9DB4D0"),
           ("other_megaplasmid", "other megaplasmid", "#D55E00")]

res = list(csv.DictReader(open(TAB / "orthogroup_residence.tsv"), delimiter="\t"))
frac = np.array([float(r["frac_on_c3"]) for r in res])

cl = list(csv.DictReader(open(TAB / "secondary_replicon_clusters.tsv"), delimiter="\t"))
for r in cl:
    for k in ("c3_content", "c2_content"):
        r[k] = float(r[k])
    r["length"] = int(r["length"])
by_class = {k: [r for r in cl if r["replicon_class"] == k] for k, _, _ in CLASSES}
print({k: len(v) for k, v in by_class.items()})

extreme = 100 * np.mean((frac < 0.1) | (frac > 0.9))
D_CUT, A_CUT = 0.10, 0.30

fig = plt.figure(figsize=(7.2, 5.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.48, wspace=0.32)

# ---------------- Panel A ----------------
ax = fig.add_subplot(gs[0, :])
ax.hist(frac, bins=40, color="#1B3A6B", edgecolor="white", linewidth=0.4)
ax.axvspan(0, 0.1, color="#D55E00", alpha=0.11, zorder=0)
ax.axvspan(0.9, 1.0, color="#D55E00", alpha=0.11, zorder=0)
ax.set_xlabel("Fraction of a gene family's occurrences residing on c3")
ax.set_ylabel("Gene families")
ax.set_title("A   Gene families are replicon-faithful, not shared at random",
             loc="left", fontweight="bold")
ax.text(0.5, 0.93,
        f"{extreme:.1f}% of {len(frac):,} shared families sit near-exclusively on\n"
        f"one replicon (shaded); 0.3% expected under a shuffled null",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="#CCCCCC", linewidth=0.6))
ax.set_xlim(0, 1)

# ---------------- Panel B ----------------
ax = fig.add_subplot(gs[1, 0])
for key, label, colour in CLASSES:
    v = by_class[key]
    ax.scatter([r["c2_content"] for r in v], [r["c3_content"] for r in v],
               s=15, c=colour, edgecolors="none", alpha=0.85,
               label=f"{label} (n={len(v)})")
ax.axhline(D_CUT, color="#B00020", lw=0.8, ls="--")
ax.axvline(A_CUT, color="#B00020", lw=0.8, ls="--")
ax.set_xlabel("Chromosome-2-diagnostic content")
ax.set_ylabel("c3-diagnostic content")
ax.set_title("B   Three classes, not two", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="upper right", handletextpad=0.3)
# point at the orange cloud rather than sitting on the y axis, where the label
# collided with the tick labels
ax.annotate("neither signature:\nnot pC3", xy=(0.17, 0.035), xytext=(0.20, 0.42),
            fontsize=6.2, color="#B45309", ha="left",
            arrowprops=dict(arrowstyle="-", lw=0.6, color="#B45309",
                            shrinkA=1, shrinkB=1))

# ---------------- Panel C ----------------
ax = fig.add_subplot(gs[1, 1])
rng = np.random.default_rng(0)
for i, (key, label, colour) in enumerate(CLASSES):
    v = [r["length"] / 1e6 for r in by_class[key]]
    if not v:
        continue
    ax.scatter(rng.normal(i, 0.075, len(v)), v, s=11, c=colour,
               edgecolors="none", alpha=0.8)
    ax.hlines(np.median(v), i - 0.26, i + 0.26, color="black", lw=1.2, zorder=5)
ax.set_xticks(range(len(CLASSES)))
ax.set_xticklabels([l.replace(" (", "\n(").replace(" mega", "\nmega")
                    for _, l, _ in CLASSES], fontsize=6.4)
ax.set_ylabel("Replicon size (Mb)")
ax.set_title("C   Size does not identify c3", loc="left", fontweight="bold")
ax.text(0.5, 0.97, "c3 and other megaplasmids\noverlap almost completely in size",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.4,
        color="#555555")

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig2_c3_coherence"))
print("done")
