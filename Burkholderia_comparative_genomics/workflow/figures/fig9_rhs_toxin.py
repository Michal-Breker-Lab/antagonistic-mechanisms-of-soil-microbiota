#!/usr/bin/env python3
"""Figure 9 - the MF6 pC3 RHS toxin across 771 Burkholderia genomes.

Panel A: per-residue alignment depth of the 2,867-aa RHS protein against the 247
         genomes that carry a related RHS protein but NOT the warhead. Depth is
         high across the RHS repeat core and falls to exactly zero at residue
         2,646 - the C-terminal 222 aa exist in no other genome. This is the
         polymorphic-toxin signature: conserved delivery scaffold, swappable tip.
Panel B: the double dissociation. The full-length RHS protein is spread across
         the genus and shows no pC3 association; the warhead is found on seven
         genomes and every one of them is pC3-positive.
Panel C: the seven carriers. Identity to the MF6 query, and the replicon the hit
         sits on - rank 3 by size in all seven, i.e. pC3.

Inputs: tables/rhs_query_coverage_profile.tsv, rhs_pC3_association.tsv,
        rhs_search_per_genome.tsv
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

PC3 = "#1B3A6B"
WARN = "#D55E00"
GREY = "#BBBBBB"
CORE = "#7FA6D9"

CT_START = 2514          # first residue of the pasted query
# QLEN and LAST_SHARED were hard-coded (2867 / 2645) in the original. They are
# now read off the coverage profile: the rebuild re-derives both exactly, but
# hard-coding them would silently mis-draw the warhead span if the query or the
# genome set ever changed again.


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


prof = read_tsv(TAB / "rhs_query_coverage_profile.tsv")
assoc = read_tsv(TAB / "rhs_pC3_association.tsv")
per = read_tsv(TAB / "rhs_search_per_genome.tsv")

QLEN = max(int(r["query_residue"]) for r in prof)
LAST_SHARED = max(int(r["query_residue"]) for r in prof
                  if int(r["n_hsps_ct_negative"]) > 0)
# denominators for the panel-B legend and the panel-A y-label
N_POS = int(assoc[0]["n_pC3pos_total"])
N_NEG = int(assoc[0]["n_pC3neg_total"])
WARHEAD = [r for r in per
           if r["query"] == "CT354_CFFIHE_03684" and r["best_tier"] == "1"]
N_CARRIER = len(WARHEAD)
N_CTNEG = len({r["accession"] for r in per
               if r["query"] == "FULL_CFFIHE_03684" and int(r["n_hits_any"]) > 0}
              - {r["accession"] for r in WARHEAD})
print(f"query length {QLEN}; last CT-negative residue {LAST_SHARED}; "
      f"carriers {N_CARRIER}; CT-negative with a FULL hit {N_CTNEG}; "
      f"pC3+ {N_POS} / pC3- {N_NEG}")

fig = plt.figure(figsize=(7.5, 7.6))
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.85, 0.95], hspace=0.62,
                      left=0.115, right=0.975, top=0.905, bottom=0.075)

# ---------------------------------------------------------------- panel A
ax = fig.add_subplot(gs[0])
x = np.array([int(r["query_residue"]) for r in prof])
y = np.array([int(r["n_hsps_ct_negative"]) for r in prof])
ax.fill_between(x, 0, y, color=CORE, linewidth=0)
ax.plot(x, y, color=PC3, lw=0.6)
ax.axvspan(LAST_SHARED + 1, QLEN, color=WARN, alpha=0.16, linewidth=0)
ax.axvline(CT_START, color="#333333", lw=0.8, ls=":")
ax.axvline(LAST_SHARED + 1, color=WARN, lw=1.0)
ax.set_xlim(1, QLEN)
ax.set_ylim(0, max(y) * 1.30)
ax.set_xlabel(f"residue position in the {QLEN:,} aa RHS protein (CFFIHE_03684)")
ax.set_ylabel(f"aligned HSPs from the\n{N_CTNEG} CT-negative genomes")
ax.annotate("pasted query starts (2514)", xy=(CT_START, max(y) * 0.72),
            xytext=(1180, max(y) * 0.68), fontsize=5.8, color="#333333",
            va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color="#333333", lw=0.7))
ax.annotate(f"warhead: residues {LAST_SHARED+1}-{QLEN}\nzero alignments from any other genome",
            xy=(2755, max(y) * 0.10),
            xytext=(1180, max(y) * 1.14), fontsize=5.8, color=WARN,
            va="center", ha="left",
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.8,
                            connectionstyle="arc3,rad=-0.18"))
ax.text(2270, max(y) * 0.07, "conserved RHS repeat core", fontsize=5.8,
        color="#12305C", ha="center")
ax.set_title("A   The C-terminal 222 aa are absent from every other genome",
             loc="left", fontweight="bold")

# ---------------------------------------------------------------- panel B
ax = fig.add_subplot(gs[1])
labels, pos_frac, neg_frac, notes = [], [], [], []
for r in assoc:
    labels.append("full-length\nRHS protein" if "FULL" in r["feature"]
                  else "C-terminal\nwarhead")
    pos_frac.append(100 * int(r["n_pC3pos_with"]) / int(r["n_pC3pos_total"]))
    neg_frac.append(100 * int(r["n_pC3neg_with"]) / int(r["n_pC3neg_total"]))
    notes.append(f"OR {r['odds_ratio']},  Fisher p = {r['fisher_p']}")
order = [1, 0]     # warhead first
labels = [labels[i] for i in order]
pos_frac = [pos_frac[i] for i in order]
neg_frac = [neg_frac[i] for i in order]
notes = [notes[i] for i in order]
yy = np.arange(len(labels))
h = 0.32
ax.barh(yy - h / 2, pos_frac, height=h, color=PC3, label=f"pC3-positive genomes (n = {N_POS})")
ax.barh(yy + h / 2, neg_frac, height=h, color=GREY, label=f"pC3-negative genomes (n = {N_NEG})")
for i in range(len(labels)):
    ax.text(max(pos_frac[i], neg_frac[i]) + 1.6, yy[i], notes[i],
            va="center", fontsize=6.0, color="#333333")
ax.set_yticks(yy); ax.set_yticklabels(labels, fontsize=7.2)
ax.set_xlabel("% of genomes carrying the feature")
ax.set_xlim(0, 52)
ax.invert_yaxis()
ax.legend(frameon=False, fontsize=6.0, loc="lower right")
ax.set_title("B   The scaffold is replicon-agnostic; the warhead is pC3-exclusive",
             loc="left", fontweight="bold")

# ---------------------------------------------------------------- panel C
ax = fig.add_subplot(gs[2])
car = [r for r in per if r["query"] == "CT354_CFFIHE_03684" and r["best_tier"] == "1"]
car.sort(key=lambda r: -float(r["best_pident"]))
names, pids, sizes = [], [], []
for r in car:
    org = r["organism"].replace("Burkholderia ", "B. ")
    names.append(f"{org}\n{r['accession']}")
    pids.append(float(r["best_pident"]))
    sizes.append(int(r["best_contig_len"]) / 1e6)
yy = np.arange(len(names))
ax.barh(yy, pids, color=PC3, height=0.55)
for i, (p, s) in enumerate(zip(pids, sizes)):
    ax.text(p - 1.2, yy[i], f"{p:.1f}%", va="center", ha="right",
            fontsize=6.2, color="white")
    ax.text(101, yy[i], f"replicon 3  ·  {s:.2f} Mb", va="center",
            fontsize=6.0, color="#333333")
ax.axvline(60, color=WARN, lw=0.9, ls="--")
ax.text(60, len(names) - 0.25, " 60% threshold", fontsize=5.8, color=WARN, va="top")
ax.set_yticks(yy); ax.set_yticklabels(names, fontsize=5.8)
ax.set_xlim(55, 118)
ax.set_xticks([60, 70, 80, 90, 100])
ax.invert_yaxis()
ax.set_xlabel("amino-acid identity to the MF6 query (100% query coverage throughout)")
ax.set_title(f"C   All {N_CARRIER} carriers, all on the third replicon",
             loc="left", fontweight="bold")
for sp in ("right", "top"):
    ax.spines[sp].set_visible(False)

fig.suptitle("An RHS polymorphic toxin of $\\it{B.\\ sola}$ pC3: conserved delivery scaffold,\n"
             f"pC3-restricted warhead shared by {N_CARRIER} genomes",
             fontsize=10, fontweight="bold", y=0.982)

ok = figstyle.save(fig, str(FIG / "fig9_rhs_toxin"))
print("OK" if ok else "FAILED vector-text check")
