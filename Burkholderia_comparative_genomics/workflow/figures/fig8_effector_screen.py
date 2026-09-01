#!/usr/bin/env python3
"""Figure 8 - InterProScan / antiSMASH / MacSyFinder screen of the B. sola pC3.

Panel A: annotation rescue. COG classified 28.1% of the pC3 clade core; adding
         InterProScan raises it to 88.7%. The residue that still matches nothing in
         any of 17 member databases is shown explicitly rather than dropped.
Panel B: biosynthetic gene cluster inventory. Region types are near-identical
         across the six replicons; ~20-23% of each pC3 is BGC.
Panel C: secretion systems called at SYSTEM level by MacSyFinder/TXSScan. No
         complete T6SS on any pC3; some carry a partial T6SS module
         (open marker) that fails quorum at 3 of 11 mandatory genes.
Panel D: the MF6 effector cassette. All genes on the minus strand, so
         transcription runs right to left: TssF, VgrG, adaptors, the PAAR-RHS-HNH
         polymorphic toxin, and a Knr4/Smi1-family immunity protein, the whole
         island flanked by transposases.

Inputs: tables/screen_ips_coverage.tsv, screen_antismash_regions.tsv,
        screen_secretion_systems.tsv, screen_secretion_rejected.tsv
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrow  # noqa: E402
import numpy as np  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

PC3 = "#1B3A6B"
ACC = "#56B4E9"
NONE_C = "#D9D9D9"
WARN = "#D55E00"

# Membership and replicon lengths are read, not hard-coded: Set B gained MF7 and
# MF7's pC3 spans two contigs, so a literal list here would silently omit it (the
# failure D16 records).
def _setb():
    rows = list(csv.DictReader(open(Path(__file__).resolve().parent.parent
                                    / "tables" / "setB_contigs.tsv", newline=""),
                               delimiter="\t"))
    order, lens = [], {}
    for r in rows:
        if r["replicon"] != "pC3":
            continue
        if r["accession"] not in order:
            order.append(r["accession"])
        lens[r["accession"]] = lens.get(r["accession"], 0) + int(r["length"])
    return order, lens


SETB, PC3_LEN = _setb()
SHORT = {"GCF_016899425.1": "MS389", "MF6": "MF6", "MF7": "MF7"}
SHORT.update({g: g.replace(".1", "").replace(".2", "") for g in SETB if g not in SHORT})


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


cov = read_tsv(TAB / "screen_ips_coverage.tsv")
regions = read_tsv(TAB / "screen_antismash_regions.tsv")
systems = read_tsv(TAB / "screen_secretion_systems.tsv")
rejected = read_tsv(TAB / "screen_secretion_rejected.tsv")

fig = plt.figure(figsize=(7.5, 8.4))
_gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.85, 0.62], hspace=0.80,
                       left=0.13, right=0.975, top=0.918, bottom=0.075)
gs_top = _gs[0].subgridspec(1, 2, wspace=0.45)

# ---------------------------------------------------------------- panel A
ax = fig.add_subplot(gs_top[0, 0])
groups = [("core", "pC3 clade core"), ("accessory", "pC3 clade accessory")]
labels, cogv, iprv, nonev = [], [], [], []
for key, lab in groups:
    rs = [r for r in cov if r["gene_set"] == key]
    n = sum(int(r["n_families"]) for r in rs)
    n_cog = sum(int(r["n_families"]) for r in rs if r["cog_status"] == "COG")
    n_ipr_extra = sum(int(r["n_interpro_entry"]) for r in rs if r["cog_status"] == "noCOG")
    labels.append(f"{lab}\n(n = {n:,})")
    cogv.append(100 * n_cog / n)
    iprv.append(100 * n_ipr_extra / n)
    nonev.append(100 - 100 * n_cog / n - 100 * n_ipr_extra / n)
y = np.arange(len(labels))
ax.barh(y, cogv, color=PC3, height=0.5, label="COG category (as in §6.4)")
ax.barh(y, iprv, left=cogv, color=ACC, height=0.5, label="InterPro entry added here")
ax.barh(y, nonev, left=np.array(cogv) + np.array(iprv), color=NONE_C, height=0.5,
        label="still unassigned")
for i in range(len(labels)):
    ax.text(cogv[i] / 2, y[i], f"{cogv[i]:.0f}", ha="center", va="center",
            fontsize=6.5, color="white")
    ax.text(cogv[i] + iprv[i] / 2, y[i], f"+{iprv[i]:.0f}", ha="center", va="center",
            fontsize=6.5, color="#123")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.8)
ax.set_xlabel("% of gene families"); ax.set_xlim(0, 100)
ax.set_title("A   InterProScan closes the annotation gap", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=5.4, ncol=3, loc="upper center",
          bbox_to_anchor=(0.5, -0.30), columnspacing=1.0, handlelength=1.2)

# ---------------------------------------------------------------- panel B
ax = fig.add_subplot(gs_top[0, 1])
ptypes = []
for r in regions:
    for p in r["products"].split(";"):
        if p and p not in ptypes:
            ptypes.append(p)
ptypes.sort()
grid = np.zeros((len(ptypes), len(SETB)))
for j, g in enumerate(SETB):
    for r in [x for x in regions if x["genome"] == g]:
        for p in r["products"].split(";"):
            if p:
                grid[ptypes.index(p), j] += 1
ax.imshow(grid, aspect="auto", cmap="Blues", vmin=0, vmax=max(3, grid.max()))
for i in range(len(ptypes)):
    for j in range(len(SETB)):
        if grid[i, j]:
            ax.text(j, i, int(grid[i, j]), ha="center", va="center", fontsize=6,
                    color="white" if grid[i, j] >= 2 else "#123")
ax.set_xticks(range(len(SETB)))
ax.set_xticklabels([SHORT[g] for g in SETB], rotation=45, ha="right", fontsize=5.8)
ax.set_yticks(range(len(ptypes))); ax.set_yticklabels(ptypes, fontsize=6)
ax.set_title("B   Biosynthetic clusters, ~20-23% of pC3", loc="left", fontweight="bold")
ax.tick_params(length=0)

# ---------------------------------------------------------------- panel C
ax = fig.add_subplot(_gs[1])
models = ["T5cSS", "Tad", "Flagellum", "T6SSi"]
NICE = {"T5cSS": "T5cSS adhesin", "Tad": "Tad pilus",
        "Flagellum": "Flagellar export", "T6SSi": "T6SS (complete)"}
w = 0.9 / len(SETB)
for k, g in enumerate(SETB):
    vals = []
    for m in models:
        vals.append(sum(1 for r in systems if r["genome"] == g and r["model"] == m))
    ax.bar(np.arange(len(models)) + (k - (len(SETB) - 1) / 2) * w, vals, width=w,
           color=plt.cm.viridis(k / max(len(SETB) - 1, 1)), label=SHORT[g])
# mark the partial T6SS
n_partial = len({r["genome"] for r in rejected if "T6SS" in r["model"]})
ax.scatter([len(models) - 1], [0.35], marker="o", s=42, facecolors="none",
           edgecolors=WARN, linewidths=1.4, zorder=5)
ax.annotate(f"partial T6SS module\nin {n_partial}/{len(SETB)} (quorum 3/11)",
            xy=(len(models) - 1, 0.35), xytext=(len(models) - 1.75, 2.4),
            fontsize=5.8, color=WARN,
            arrowprops=dict(arrowstyle="->", color=WARN, lw=0.8))
ax.set_xticks(range(len(models)))
ax.set_xticklabels([NICE[m] for m in models], fontsize=7.0)
ax.set_ylabel("systems per pC3")
ax.set_xlim(-0.55, len(models) - 0.45)
ax.set_title("C   Secretion systems, system-level calls", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=5.8, ncol=len(SETB), loc="upper center",
          bbox_to_anchor=(0.5, -0.14), columnspacing=1.2, handlelength=1.2)

# ---------------------------------------------------------------- panel D
ax = fig.add_subplot(_gs[2])
GENES = [
    ("CFFIHE_04289", 1058182, 1058706, "Knr4/Smi1\nimmunity", "#009E73"),
    ("CFFIHE_04290", 1058716, 1062963, "PAAR-RHS-HNH\ntoxin", WARN),
    ("CFFIHE_04291", 1062985, 1063443, "", "#BBBBBB"),
    ("CFFIHE_04292", 1063509, 1063952, "", "#BBBBBB"),
    ("CFFIHE_04293", 1063996, 1064865, "Ank", "#BBBBBB"),
    ("CFFIHE_04294", 1064865, 1067864, "VgrG", PC3),
    ("CFFIHE_04295", 1068012, 1069778, "TssF", PC3),
    ("CFFIHE_04296", 1070150, 1070602, "transposase", "#8C6D31"),
    ("CFFIHE_04297", 1070640, 1070741, "", "#8C6D31"),
]
x0 = GENES[0][1] - 800
x1 = GENES[-1][2] + 800
for lt, s, e, lab, col in GENES:
    # every gene is on the minus strand, so arrows point left
    ax.add_patch(FancyArrow(e, 0, s - e, 0, width=0.30, head_width=0.5,
                            head_length=min(320, (e - s) * 0.45),
                            length_includes_head=True, color=col, linewidth=0))
    if lab:
        ax.text((s + e) / 2, 0.55, lab, ha="center", va="bottom", fontsize=5.8)
ax.annotate("", xy=(x0 + 1200, -0.85), xytext=(x0 + 4200, -0.85),
            arrowprops=dict(arrowstyle="->", lw=1.0, color="#333333"))
ax.text(x0 + 4400, -0.85, "transcription", fontsize=5.8, va="center")
ax.set_xlim(x0, x1); ax.set_ylim(-1.3, 1.5)
ax.set_yticks([])
ax.set_xlabel("position on MF6 pC3 (bp)")
ax.set_title("D   The MF6 pC3 effector cassette, flanked by transposases",
             loc="left", fontweight="bold")
for sp in ("left", "right", "top"):
    ax.spines[sp].set_visible(False)

fig.suptitle("Screening the $\\it{B.\\ sola}$ pC3: an orphan VgrG-PAAR effector cassette,\n"
             "toxin-antitoxin loci, and an antifungal biosynthetic arsenal",
             fontsize=10, fontweight="bold", y=0.985)

ok = figstyle.save(fig, str(FIG / "fig8_effector_screen"))
print("OK" if ok else "FAILED vector-text check")
