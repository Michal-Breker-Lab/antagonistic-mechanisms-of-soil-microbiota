#!/usr/bin/env python3
"""Figure 7 - functional composition of the B. sola clade pC3 core.

Panel A: COG functional-category composition of the pC3 clade core (833 families,
         5 genomes) beside the chromosome-1 clade core (2,847 families). Percentages
         are of COG-ANNOTATED families; a family carrying several categories is
         counted in each, so columns sum to slightly over 100%.
Panel B: descriptive effect sizes for the same contrast - log2 odds ratio with 95%
         CI (Haldane-Anscombe corrected). These are NOT inferential: the two sets
         are two replicons of the same five genomes, not independent samples, so
         the Fisher null of random gene-to-replicon assignment is not a real null.
         BH-corrected q values are marked only to show which contrasts are robust.
Panel C: annotation coverage, reported as a result rather than silently dropped.
         The dominant feature of the pC3 core is that most of it carries a named
         product but no COG category at all.
Panel D: KEGG modules whose completeness rises when pC3 is added to the rest of
         the MF6 genome - i.e. what pC3 contributes metabolically. Stepwise
         completeness, anvi'o definition, from Bakta KO calls (see TOOLS.md).

Inputs: tables/setB_cog_profile.tsv, setB_cog_enrichment.tsv,
        setB_product_buckets.tsv, setB_pC3_module_contribution.tsv
"""
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

PC3 = "#1B3A6B"      # pC3, matching C3_COLORS[True] elsewhere in the report
CHR1 = "#8C6D31"     # chromosome 1, matching fig4
UP, DOWN = "#0072B2", "#D55E00"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


prof = read_tsv(TAB / "setB_cog_profile.tsv")
enr = read_tsv(TAB / "setB_cog_enrichment.tsv")
buck = read_tsv(TAB / "setB_product_buckets.tsv")
modc = read_tsv(TAB / "setB_pC3_module_contribution.tsv")

# ---- contrast A, COG-annotated families only -------------------------------
rowsA = [r for r in enr
         if r["mode"] == "COG_annotated_only" and r["contrast"] == "A_pC3core_vs_chr1core"
         and r["cog_category"] != "-"]
# drop categories that are essentially absent from both sets
rowsA = [r for r in rowsA if int(r["n_a"]) + int(r["n_b"]) >= 5]
rowsA.sort(key=lambda r: float(r["log2_OR"]))

cats = [r["cog_category"] for r in rowsA]
names = [r["cog_name"] for r in rowsA]
pct_a = np.array([float(r["pct_a"]) for r in rowsA])
pct_b = np.array([float(r["pct_b"]) for r in rowsA])
lor = np.array([float(r["log2_OR"]) for r in rowsA])
lo = np.array([float(r["ci_low"]) for r in rowsA])
hi = np.array([float(r["ci_high"]) for r in rowsA])
q = np.array([float(r["q_BH"]) for r in rowsA])

# 7.5 x 8.4 keeps the whole figure plus its caption inside one LaTeX text block
# (6.5 in x 9 in). Taller aspects overflow the float and pandoc then distorts the
# figure to fit, which silently squashes every label.
fig = plt.figure(figsize=(7.5, 8.4))
_gs = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.32,
                       left=0.20, right=0.975, top=0.925, bottom=0.075)
# the bottom row needs a much wider gutter: panel D's module names are long and
# would otherwise be drawn straight over panel C's bars
gs_top = _gs[0].subgridspec(1, 2, wspace=0.42)
gs_bot = _gs[1].subgridspec(1, 2, wspace=1.05)

# ---------------------------------------------------------------- panel A
ax = fig.add_subplot(gs_top[0, 0])
y = np.arange(len(cats))
h = 0.38
# Set sizes come from the enrichment table's own totals, not from the caption:
# the rebuild changed them (2,847 -> 2,833) and a hard-coded legend would lie.
N_A = int(rowsA[0]["total_a"])
N_B = int(rowsA[0]["total_b"])
ax.barh(y + h / 2, pct_a, height=h, color=PC3, label=f"pC3 clade core (n = {N_A:,})")
ax.barh(y - h / 2, pct_b, height=h, color=CHR1,
        label=f"chromosome 1 clade core (n = {N_B:,})")
ax.set_yticks(y)
ax.set_yticklabels([f"{c}  {n}" for c, n in zip(cats, names)], fontsize=6.2)
ax.set_xlabel("% of COG-annotated families")
ax.set_title("A   Functional composition", loc="left", fontweight="bold")
# categories are sorted by effect size, so the two top rows (W, B) have short bars
# on both sides - the only reliably empty corner of this panel
ax.legend(frameon=True, facecolor="white", edgecolor="none", framealpha=0.9,
          fontsize=6.4, loc="upper right", bbox_to_anchor=(1.0, 1.0))
ax.set_ylim(-0.7, len(cats) - 0.3)

# ---------------------------------------------------------------- panel B
ax = fig.add_subplot(gs_top[0, 1])
colors = [UP if v > 0 else DOWN for v in lor]
ax.hlines(y, lo, hi, color="#999999", linewidth=0.9, zorder=1)
ax.scatter(lor, y, s=16, c=colors, zorder=3, edgecolors="none")
ax.axvline(0, color="#333333", linewidth=0.8, linestyle="--", zorder=0)
for i, qq in enumerate(q):
    if qq < 0.05:
        ax.text(hi[i] + 0.25, y[i], "*", va="center", fontsize=9, fontweight="bold")
ax.set_yticks(y)
ax.set_yticklabels(cats, fontsize=6.5)
ax.set_xlabel("log$_2$ odds ratio (pC3 core vs chr1 core)")
ax.set_title("B   Effect size, descriptive", loc="left", fontweight="bold")
ax.set_ylim(-0.7, len(cats) - 0.3)
ax.text(0.03, 0.985, "* BH q < 0.05\nright = enriched on pC3",
        transform=ax.transAxes, fontsize=6.2, va="top", color="#444444")

# ---------------------------------------------------------------- panel C
ax = fig.add_subplot(gs_bot[0, 0])
SETS = [("chr1_clade_core", "chr1 clade core"),
        ("pC3_clade_accessory", "pC3 clade accessory"),
        ("pC3_clade_core", "pC3 clade core")]
# Buckets are computed here (not by subtracting the coverage table) so that they
# are mutually exclusive: a COG-assigned family whose product happens to be a DUF
# must not be counted twice.
fams = read_tsv(TAB / "setB_functional_families.tsv")
import re as _re


def _bucket(r):
    if r["cog_categories"].strip():
        return "COG category assigned"
    p = r["product"].lower()
    if "hypothetical" in p:
        return "hypothetical protein"
    if _re.search(r"\bduf\d+", p) or "uncharacteri" in p or "domain-containing protein" in p:
        return "domain-containing / DUF only"
    return "named product, no COG"


counts = {}
for r in fams:
    key = r["gene_set"].split("/")[0]
    counts.setdefault(key, []).append(_bucket(r))
cov = {k: {"n_families": len(v)} for k, v in counts.items()}

BUCKETS = [("COG category assigned", "#1B3A6B"),
           ("named product, no COG", "#56B4E9"),
           ("domain-containing / DUF only", "#F0E442"),
           ("hypothetical protein", "#D9D9D9")]
yy = np.arange(len(SETS))
for i, (key, lab) in enumerate(SETS):
    n = len(counts[key])
    from collections import Counter as _C
    cc = _C(counts[key])
    vals = [100 * cc.get(b, 0) / n for b, _ in BUCKETS]
    left = 0.0
    for val, (lab2, col) in zip(vals, BUCKETS):
        ax.barh(yy[i], val, left=left, color=col, height=0.6,
                label=lab2 if i == 0 else None, edgecolor="white", linewidth=0.4)
        if val > 7:
            ax.text(left + val / 2, yy[i], f"{val:.0f}", ha="center", va="center",
                    fontsize=6.2, color="white" if col == "#1B3A6B" else "#222222")
        left += val
ax.set_yticks(yy)
ax.set_yticklabels([f"{l}\n(n = {cov[k]['n_families']:,})" for k, l in SETS],
                   fontsize=6.8)
ax.set_xlabel("% of gene families")
ax.set_xlim(0, 100)
ax.set_title("C   How much is actually annotated", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.0, ncol=2, loc="upper center",
          bbox_to_anchor=(0.5, -0.28))

# ---------------------------------------------------------------- panel D
ax = fig.add_subplot(gs_bot[0, 1])
mc = [r for r in modc if r["pC3_raises_completeness"] == "1"]
for r in mc:
    r["_gain"] = float(r["whole_genome"]) - float(r["rest_of_genome"])
mc.sort(key=lambda r: r["_gain"], reverse=True)
mc = mc[:12][::-1]
yy = np.arange(len(mc))
rest = np.array([float(r["rest_of_genome"]) for r in mc]) * 100
gain = np.array([r["_gain"] for r in mc]) * 100
ax.barh(yy, rest, color="#CCCCCC", height=0.62, label="rest of the MF6 genome")
ax.barh(yy, gain, left=rest, color=PC3, height=0.62, label="added by pC3")


def short(nm):
    nm = nm.split(",")[0]
    return nm[:38] + ("…" if len(nm) > 38 else "")


ax.set_yticks(yy)
ax.set_yticklabels([f"{short(r['name'])}" for r in mc], fontsize=6.0)
ax.set_xlabel("KEGG module stepwise completeness (%)")
ax.set_xlim(0, 105)
ax.set_title("D   What pC3 adds metabolically", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.2, loc="lower right")

fig.suptitle("Functional composition of the $\\it{B.\\ sola}$ pC3 core: a regulatory,\n"
             "transport and peripheral-catabolic replicon, mostly outside COG",
             fontsize=10, fontweight="bold", y=0.985)

ok = figstyle.save(fig, str(FIG / "fig7_functional_composition"))
print("OK" if ok else "FAILED vector-text check")
