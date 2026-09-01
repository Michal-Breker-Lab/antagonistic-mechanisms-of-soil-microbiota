#!/usr/bin/env python3
"""Figure 4 - the c3 pangenome.

Panel A: rarefaction. Pan-genome and core-genome size as genomes are added,
         with Heaps'-law fit. gamma > 0 means an OPEN pangenome -- new genomes
         keep contributing new genes.
Panel B: gene-frequency spectrum (core / soft-core / shell / cloud).
Panel C: clean-mode sensitivity. Panaroo's strict mode discards suspected
         contaminants present in <5% of genomes, moderate uses 1%. For a
         megaplasmid pangenome that difference is not cosmetic, so all three
         modes are reported side by side.

Input: Panaroo gene_presence_absence.Rtab (genes x genomes, binary).
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
rng = np.random.default_rng(0)

MODES = [("c3_moderate", "moderate (primary)", "#1B3A6B"),
         ("c3_sensitive", "sensitive", "#56B4E9"),
         ("chr1_strict", "chromosome 1 (strict)", "#8C6D31")]


def load_rtab(path):
    with open(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        genomes = hdr[1:]
        rows = [np.array([int(x) for x in r[1:]], dtype=bool) for r in rd if len(r) > 1]
    return np.array(rows), genomes


def rarefy(M, n_perm=25):
    """-> (pan_mean, core_mean) as a function of number of genomes."""
    n_g = M.shape[1]
    pan = np.zeros((n_perm, n_g)); core = np.zeros((n_perm, n_g))
    for p in range(n_perm):
        order = rng.permutation(n_g)
        seen = np.zeros(M.shape[0], dtype=bool)
        cnt = np.zeros(M.shape[0], dtype=int)
        for i, g in enumerate(order):
            col = M[:, g]
            seen |= col
            cnt += col
            pan[p, i] = seen.sum()
            core[p, i] = (cnt == (i + 1)).sum()
    return pan.mean(0), core.mean(0)


avail = {}
for key, label, colour in MODES:
    p = TAB / f"{key}_gene_presence_absence.Rtab"
    if p.exists():
        avail[key] = load_rtab(p)
        print(f"{key}: {avail[key][0].shape[0]} genes x {len(avail[key][1])} genomes")
    else:
        print(f"{key}: MISSING ({p.name})")

if "c3_moderate" not in avail:
    sys.exit("c3_moderate pangenome not available - run Stage 7 first")

fig = plt.figure(figsize=(7.2, 5.4))
gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.30)

# ---------------- Panel A: rarefaction ----------------
ax = fig.add_subplot(gs[0, :])
M, genomes = avail["c3_moderate"]
pan, core = rarefy(M)
x = np.arange(1, len(pan) + 1)
ax.plot(x, pan, color="#1B3A6B", lw=1.6, label="pan-genome")
ax.plot(x, core, color="#D55E00", lw=1.6, label="core genome")

# Heaps' law on the pan curve: n_genes ~ k * N^gamma.
# np.polyfit returns [slope, intercept]; the SLOPE in log-log space is gamma.
mask = x >= 3
gamma_exp, log_k = np.polyfit(np.log(x[mask]), np.log(pan[mask]), 1)
ax.plot(x, np.exp(log_k) * x ** gamma_exp, color="#999999", lw=0.9, ls="--",
        label=f"Heaps' fit, $\\gamma$={gamma_exp:.2f}")
ax.set_xlabel("Genomes sampled")
ax.set_ylabel("Gene families")
ax.set_title(f"A   c3 pangenome rarefaction (n = {len(genomes)} genomes)",
             loc="left", fontweight="bold")
ax.legend(frameon=False, loc="center right")
ax.text(0.98, 0.34,
        ("open pangenome: new genomes keep adding genes"
         if gamma_exp > 0.05 else "closed pangenome"),
        transform=ax.transAxes, ha="right", fontsize=7, color="#555555")

# ---------------- Panel B: frequency spectrum ----------------
ax = fig.add_subplot(gs[1, 0])
freq = M.sum(1) / M.shape[1]
ax.hist(freq, bins=40, color="#1B3A6B", edgecolor="white", linewidth=0.4)
ax.set_xlabel("Fraction of genomes carrying the gene family")
ax.set_ylabel("Gene families")
ax.set_title("B   c3 gene-frequency spectrum", loc="left", fontweight="bold")
cats = [("core >=95%", (freq >= 0.95).sum()),
        ("shell 15-95%", ((freq >= 0.15) & (freq < 0.95)).sum()),
        ("cloud <15%", (freq < 0.15).sum())]
ax.text(0.52, 0.95, "\n".join(f"{n:,}  {c}" for c, n in cats),
        transform=ax.transAxes, va="top", fontsize=7,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                  edgecolor="#CCCCCC", linewidth=0.6))

# ---------------- Panel C: clean-mode sensitivity ----------------
ax = fig.add_subplot(gs[1, 1])
labels, cores, totals, cols = [], [], [], []
for key, label, colour in MODES:
    if key not in avail:
        continue
    Mm, gg = avail[key]
    f = Mm.sum(1) / Mm.shape[1]
    labels.append(label.replace(" (", "\n("))
    cores.append((f >= 0.95).sum())
    totals.append(Mm.shape[0])
    cols.append(colour)
xx = np.arange(len(labels))
# log scale: core (tens) and total (tens of thousands) differ by ~3 orders of
# magnitude, so a linear axis renders the core bars invisible
ax.bar(xx - 0.16, totals, width=0.32, color="#D9D9D9", edgecolor="white",
       linewidth=0.5, label="total families")
ax.bar(xx + 0.16, cores, width=0.32, color=cols, edgecolor="white",
       linewidth=0.5, label="core (>=95%)")
ax.set_yscale("log")
ax.set_xticks(xx); ax.set_xticklabels(labels, fontsize=6.4)
ax.set_ylabel("Gene families (log scale)")
ax.set_title("C   Clean-mode sensitivity", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6.4, loc="upper left")
for i, (t, c) in enumerate(zip(totals, cores)):
    ax.text(i - 0.16, t * 1.15, f"{t:,}", ha="center", va="bottom", fontsize=5.8)
    ax.text(i + 0.16, c * 1.15, f"{c:,}", ha="center", va="bottom", fontsize=5.8)
ax.set_ylim(top=max(totals) * 6)

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig4_c3_pangenome"))
print("done")
