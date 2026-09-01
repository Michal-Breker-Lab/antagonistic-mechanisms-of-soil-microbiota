#!/usr/bin/env python3
"""Figure 11 - where the 10 candidate toxin / immunity proteins sit on the
chromosome-1 phylogeny.

Same chromosome-1 tree as Figure 4, which is left untouched; this is a separate
figure because Figure 1's job is the pC3 story and ten more rings would bury it.

Panel A: the tree with one ring per locus, ordered inner -> outer by how patchy
         the locus is, so the informative ones sit nearest the tree. Every locus
         gets its own colour. Two rings sit between the tree and the toxin block:
         pC3 presence, then isolation source.

Revised 2026-08-13 at Moshe's request:
  * For the two RHS proteins searched with BOTH a full-length query and a
    C-terminal tip query, the ring now shows the TIP only. The tip is the
    polymorphic warhead; the full-length hit is dominated by the conserved RHS
    scaffold, so a scaffold hit does not establish that the same warhead is
    present. For MF6_004284 this matters -- 8 carriers by the scaffold, 4 by the
    tip. Panel B was trimmed the same way on 2026-08-13: the two full-length
    scaffold rows are gone, so the whole figure now speaks about warheads only
    and no panel can be read as claiming a scaffold hit is a toxin match. The
    scaffold numbers remain in toxin12_carriers.tsv and toxin12_phylo.tsv.
  * Every locus now gets a distinct hue. Previously the four near-universal loci
    (>=400 carriers) were drawn grey to signal that an almost-solid ring carries
    little distributional information. That signal now lives only in the carrier
    count printed in the key.
  * The isolation-source ring from Figure 1 is back, using the same colours and
    the same host_categories.tsv, so the two figures can be read against
    each other.
Panel B: the two RHS cassettes side by side across the union of their carriers,
         each as its C-terminal warhead plus the dark protein next to it. This is
         the result the rings cannot show - MF6 has TWO RHS polymorphic toxins on
         pC3, and they have different carrier sets.
Panel C: phylogenetic clustering per locus, against the pC3-positive null.
         Negative z = carriers closer together than chance. Loci with more
         carriers than there are pC3-positive tips have no such null and are
         shown as "n/a" rather than silently omitted.

Inputs: tables/chr1_core.treefile, c3_calls_all_genomes.tsv,
        toxin12_carriers.tsv, toxin12_locus_summary.tsv, toxin12_phylo.tsv,
        host_categories.tsv
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import radialtree as rt  # noqa: E402
import tree_display as td  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)

# ONE question, one encoding: which protein is in which genome. Each ring is
# filled in its own protein colour where the genome carries that locus, and the
# same grey where it does not. An earlier version also coloured fills by the
# replicon the hit sat on; that answered a second question at the cost of making
# the first one hard to read, and was removed. Replicon placement is still
# recorded, in tables/toxin12_carrier_replicon.tsv.
# Ten distinct hues (Paul Tol muted, colourblind-safe) -- one per locus, no
# recycling and no greying, so ring colour identifies the protein uniquely.
PROT_COLOR = ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
              "#DDCC77", "#CC6677", "#882255", "#AA4499", "#EE7733"]
# still used to decide which rings get a widened wedge, not to decide colour
MAX_HITS_COLOURED = 400

# The rings to draw, one per protein. Where a locus was searched with both a
# full-length query and a C-terminal tip query, the TIP is the ring -- see the
# module docstring. Panel B keeps both.
RING_LOCI = ["MF6_001079", "MF6_002734", "MF6_003184",
             "MF6_003684_Cterm", "MF6_003686", "MF6_003843",
             "MF6_004284_Cterm", "MF6_004285", "MF6_004947", "MF6_006318"]


def base_locus(locus):
    """Ring key -> the row that carries this locus's metadata in
    toxin12_locus_summary.tsv, which is keyed by the full-length locus."""
    return locus[:-len("_Cterm")] if locus.endswith("_Cterm") else locus


def ring_label(locus):
    return (base_locus(locus) + "  C-term tip") if locus.endswith("_Cterm") \
        else locus
PRESENT = "#1B3A6B"      # panel B, and the inner pC3 ring
ABSENT = "#EDEDED"       # absent -- the SAME grey everywhere, no second shade
CLUSTER = "#009E73"      # panel C only: clustered on the tree
NOTSIG = "#B0B6BC"       # panel C only: not significant
BRANCH = "#9AA0A6"       # tree branches
STAR = "#D55E00"         # MF6, and nothing else


def draw_ring_minwidth(ax, lay, values, colors, r0, r1, min_deg=0.0,
                       default=ABSENT):
    """One wedge per tip, but with a floor on how narrow a PRESENT wedge may be.

    Why this exists: 763 tips over 352 degrees is 0.46 deg per tip, so a locus
    carried by six genomes covers 6.9 deg -- about 2% of the circle, drawn two
    or three pixels wide. The restricted loci, which are the interesting ones,
    were effectively invisible while the four near-universal loci dominated the
    figure. Widening present wedges costs angular precision, so it is applied
    only to the sparse loci, where there is nothing adjacent to overwrite.
    """
    from matplotlib.patches import Wedge as _W
    tips = lay.tips
    angs = sorted(lay.tip_angle(t) for t in tips)
    step = (angs[-1] - angs[0]) / max(len(tips) - 1, 1)
    # single background annulus, then paint the present tips on top
    ax.add_patch(_W((0, 0), r1, np.rad2deg(angs[0] - step / 2),
                    np.rad2deg(angs[-1] + step / 2), width=r1 - r0,
                    facecolor=default, edgecolor="none", zorder=3))
    half = max(np.rad2deg(step) / 2, min_deg / 2)
    for t in tips:
        v = values.get(t.name)
        if not v:
            continue
        a = np.rad2deg(lay.tip_angle(t))
        ax.add_patch(_W((0, 0), r1, a - half, a + half, width=r1 - r0,
                        facecolor=colors.get(v, default), edgecolor="none",
                        zorder=4))


# ------------------------------------------------------------------ data
loci = {}
for r in csv.DictReader(open(TAB / "toxin12_locus_summary.tsv"), delimiter="\t"):
    loci[r["locus"]] = r

phylo = {}
for r in csv.DictReader(open(TAB / "toxin12_phylo.tsv"), delimiter="\t"):
    phylo[r["locus"]] = r

carriers = defaultdict(set)
org = {}
_searched = set()
for r in csv.DictReader(open(TAB / "toxin12_carriers.tsv"), delimiter="\t"):
    if r["tier1_carrier"] == "True":
        carriers[r["locus"]].add(r["accession"])
    org[r["accession"]] = r["organism"]
    _searched.add(r["accession"])
# The "n/N" denominator is the number of genomes SEARCHED, read from the table
# rather than written in. It was hard-coded as 772 (the original run's 771 NCBI
# + MF6); the rebuild searched 773, so every ring caption was understating the
# denominator by one genome.
N_GENOMES = len(_searched)
print(f"toxin loci searched across {N_GENOMES} genomes")

# which replicon does each carrier's best hit sit on, in that genome?
carrier_rank = defaultdict(dict)
for r in csv.DictReader(open(TAB / "toxin12_carrier_replicon.tsv"), delimiter="\t"):
    carrier_rank[r["locus"]][r["accession"]] = min(int(r["best_contig_rank"]), 4)

c3 = {}
for r in csv.DictReader(open(TAB / "c3_calls_all_genomes.tsv"), delimiter="\t"):
    c3[r["accession"]] = (r["c3_present"] == "True")
c3["MF6"] = True

# isolation source -- the SAME file and the same controlled vocabulary Figure 1
# uses, so the two figures are directly comparable. MF6 is a lab isolate with no
# BioSample and is carried there as a curated row.
host = {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    host[r["accession"]] = r["host_category"]
if "MF6" not in host:
    print("WARNING: MF6 missing from host_categories.tsv - draws as unknown")

tree = Phylo.read(str(TAB / "chr1_core.treefile"), "newick")
# Rooted on the ingroup stem, exactly as Figure 4, via the shared conventions in
# tree_display -- the two figures draw the same tree and must not disagree about
# how it is oriented. Rooting is a DISPLAY convention only; no statistic in this
# figure depends on it (the clan test and the patristic distances are both
# rooting-independent).
#
# The near-clone COLLAPSING used in Figures 4 and 10 is deliberately NOT applied
# here. There, every collapsed wedge is uniform for the trait being drawn (a
# carrier is never collapsed), so a wedge's ring cell is a real measured value.
# Here there are ten locus rings and they are NOT uniform inside those blocks --
# three of the ten sit in over 80% of all genomes -- so collapsing would replace
# genuine per-genome variation with a "mixed" cell in seven rings out of ten.
# The rooting is free; the collapsing would cost information.
tree.root_at_midpoint()          # Bio.Phylo needs a root to walk from
print("display conventions:")
ingroup, outgroup = td.ingroup_root(tree, c3, org)
tree.ladderize()
tips = tree.get_terminals()
_n_pc3_tips = sum(1 for t in tips if c3.get(t.name))
print(f"pC3-positive tips on the tree: {_n_pc3_tips}")

# rings ordered by how patchy the locus is: fewest carriers innermost. Counted
# from the carrier sets directly rather than from n_carriers_771, because the
# _Cterm rings have their own carrier sets that the locus summary does not hold
# (MF6_004284: 8 by the scaffold, 4 by the tip).
ORDER = sorted(RING_LOCI, key=lambda k: len(carriers[k]))
for k in ORDER:
    print(f"  ring {k:22s} {len(carriers[k]):4d}/{N_GENOMES} carriers")

# ------------------------------------------------------------------ figure
# taller than the pre-2026-08-13 version: the isolation-source legend needs a
# full-width strip at the bottom, as in Figure 1
fig = plt.figure(figsize=(7.5, 9.8))
# bottom is deep: panel B's rotated genome labels are long, and the isolation
# -source legend sits below them
gs = fig.add_gridspec(2, 2, height_ratios=[1.42, 1.0], width_ratios=[1.62, 1.0],
                      hspace=0.13, wspace=0.42,
                      left=0.045, right=0.975, top=0.915, bottom=0.215)

# ---------------------------------------------------------------- panel A
axA = fig.add_subplot(gs[0, 0])
axA.set_aspect("equal"); axA.axis("off")
lay = rt.RadialLayout(tree, start_deg=92, extent_deg=352)
rt.draw_tree(axA, lay, lw=0.3, color=BRANCH)

R0 = 1.04
_ext = [[rt.polar_to_xy(lay.tip_radius(t), lay.tip_angle(t)),
         rt.polar_to_xy(R0, lay.tip_angle(t))] for t in tips]
axA.add_collection(LineCollection(_ext, colors="#CCCCCC", linewidths=0.2,
                                  linestyles=(0, (1, 2)), zorder=1))

# two genome-property rings first, then the toxin block: pC3 presence anchors
# the toxin rings against the replicon, and isolation source sits next to it so
# the habitat signal is read against the tree rather than across ten rings
rt.draw_ring(axA, lay, {t.name: c3.get(t.name, False) for t in tips},
             {True: PRESENT, False: ABSENT}, R0, R0 + 0.038, default=ABSENT)
_RHOST = R0 + 0.038 + 0.010
rt.draw_ring(axA, lay, {t.name: host.get(t.name, "unknown") for t in tips},
             figstyle.HOST_COLORS, _RHOST, _RHOST + 0.038,
             default=figstyle.HOST_COLORS["unknown"])
_n_host = sum(1 for t in tips
              if host.get(t.name, "unknown") not in ("unknown", None))
print(f"tips with a known isolation source: {_n_host}/{len(tips)}")

W, GAPR = 0.030, 0.008
ring_r, prot_col = {}, {}
r = _RHOST + 0.038 + 0.022
for _pi, locus in enumerate(ORDER):
    n_tot = len(carriers[locus])
    # one hue per locus, no recycling: ring colour identifies the protein
    prot_col[locus] = PROT_COLOR[_pi % len(PROT_COLOR)]
    # sparse loci get a minimum wedge width so a single carrier is visible;
    # the near-universal ones are left at true width, where widening would
    # merge genuinely separate blocks
    draw_ring_minwidth(axA, lay,
                       {t.name: (t.name in carriers[locus]) for t in tips},
                       {True: prot_col[locus]}, r, r + W,
                       min_deg=2.6 if n_tot < MAX_HITS_COLOURED else 0.0,
                       default=ABSENT)
    ring_r[locus] = (r, r + W)
    r += W + GAPR

# Ticks marking every genome that carries at least one of the five RESTRICTED
# loci (<=10 carriers). Without these the eye goes straight to the four
# near-solid outer rings, which are the uninformative ones -- the loci that
# matter cover 2% of the circumference and read as empty.
RESTRICTED = [k for k in ORDER if len(carriers[k]) <= 10]
_any_restricted = set().union(*(carriers[k] for k in RESTRICTED))
_marked = 0
for t in tips:
    if t.name in _any_restricted:
        a = lay.tip_angle(t)
        x0, y0 = rt.polar_to_xy(r + 0.012, a)
        x1, y1 = rt.polar_to_xy(r + 0.040, a)
        axA.plot([x0, x1], [y0, y1], color=STAR, lw=1.1, solid_capstyle="butt",
                 zorder=8)
        _marked += 1
print(f"restricted loci: {RESTRICTED}")
print(f"tips carrying >=1 restricted locus, ticked: {_marked}")

# mark MF6
for t in tips:
    if t.name == "MF6":
        a = lay.tip_angle(t)
        x, y = rt.polar_to_xy(r + 0.045, a)
        axA.scatter([x], [y], s=70, marker="*", c=STAR, edgecolors="black",
                    linewidths=0.5, zorder=9)
        xt, yt = rt.polar_to_xy(r + 0.10, a)
        axA.text(xt, yt, "MF6", fontsize=7.5, fontweight="bold", color=STAR,
                 ha="left" if np.cos(a) > 0 else "right", va="center", zorder=9)

lim = r + 0.17
axA.set_xlim(-lim, lim); axA.set_ylim(-lim, lim)

axA.set_title(f"A   Ten candidate toxin / immunity loci on the {len(tips)}-tip "
              "chromosome-1 tree,\n      with pC3 presence and isolation source",
              loc="left", fontweight="bold", fontsize=9.0)

# Ring key gets its own axes -- drawn inside panel A it overlapped the tree.
# The vertical budget is tight (11 ring entries + a 5-row replicon legend + two
# notes must fit in 0..1), so the spacings below are chosen to total < 1.0; an
# earlier version overflowed and collided with panel C.
axK = fig.add_subplot(gs[0, 1])
axK.axis("off"); axK.set_xlim(0, 1); axK.set_ylim(0, 1)
ky = 0.995
axK.text(0.0, ky, "rings, inner \u2192 outer", fontsize=7.2, fontweight="bold",
         va="top")
ky -= 0.040
axK.add_patch(plt.Rectangle((0.0, ky - 0.010), 0.070, 0.020,
                            facecolor=PRESENT, lw=0))
axK.text(0.095, ky, "pC3 present", fontsize=6.3, va="center")
ky -= 0.028
axK.add_patch(plt.Rectangle((0.0, ky - 0.010), 0.070, 0.020,
                            facecolor="#BBBBBB", lw=0))
axK.text(0.095, ky, "isolation source (key below)", fontsize=6.3, va="center")
ky -= 0.046
for locus in ORDER:
    info = loci[base_locus(locus)]
    low = info["low_information"] == "True"
    axK.add_patch(plt.Rectangle((0.0, ky - 0.010), 0.070, 0.020,
                                facecolor=prot_col[locus], lw=0))
    n_tot = len(carriers[locus])                     # includes MF6 itself
    axK.text(0.095, ky + 0.010,
             f"{ring_label(locus)}   {info['replicon']}   {n_tot}/{N_GENOMES}",
             fontsize=6.0, va="center", fontweight="bold" if not low else "normal")
    prod = info["pgap_product"]
    prod = prod if len(prod) <= 33 else prod[:32] + "\u2026"
    axK.text(0.095, ky - 0.011, prod, fontsize=5.2, va="center", color="#555555")
    ky -= 0.050

ky -= 0.006
axK.add_patch(plt.Rectangle((0.0, ky - 0.010), 0.070, 0.020,
                            facecolor=ABSENT, lw=0))
axK.text(0.095, ky, "absent (all rings)", fontsize=6.3, va="center")

ky -= 0.040
axK.plot([0.018, 0.018], [ky - 0.009, ky + 0.009], color=STAR, lw=1.4,
         solid_capstyle="butt")
axK.text(0.095, ky, f"tick on the tree: one of the {_marked} tips carrying",
         fontsize=5.7, color=STAR, va="center")
ky -= 0.019
axK.text(0.095, ky, f"\u22651 of the five restricted loci ({len(_any_restricted)} genomes in all;",
         fontsize=5.7, color=STAR, va="center")
ky -= 0.019
# NOT "dereplicated away" any more -- the rebuild does not dereplicate (D11).
# The gap between carriers and ticked tips is MF7, which D13 removed from the
# chromosome-1 tree for carrying 93.96% gaps in the core alignment.
axK.text(0.095, ky, f"{len(_any_restricted) - _marked} not on the tree, D13)",
         fontsize=5.7, color=STAR, va="center")

ky -= 0.040
axK.text(0.0, ky, f"read the n/{N_GENOMES} count, not the amount of colour: an\n"
         "almost-solid ring is near-universal and says little.\n"
         "For the two RHS loci both the ring and panel B are the\n"
         "C-terminal TIP — the polymorphic warhead, not the\n"
         "conserved scaffold. Scaffold hits: toxin12_carriers.tsv.\n"
         "The replicon after each name is where that gene sits\n"
         "IN MF6, not a claim about the marked genome \u2014 see\n"
         "toxin12_carrier_replicon.tsv. Locus tags are PGAP.",
         fontsize=5.3, color="#666666", va="top", linespacing=1.45)

# ---------------------------------------------------------------- panel B
axB = fig.add_subplot(gs[1, 0])
# Trimmed 2026-08-13 to match panel A: tips only, no full-length scaffold rows.
# Dropping them also shrinks the column set, because the union of carriers is
# taken over these loci and the scaffolds pulled in genomes that carry the
# conserved RHS repeat without the matching warhead.
CASS = [("MF6_003684_Cterm", "MF6_003684  C-term tip"),
        ("MF6_003686", "MF6_003686  dark, adjacent"),
        ("MF6_004284_Cterm", "MF6_004284  C-term tip"),
        ("MF6_004285", "MF6_004285  dark, adjacent")]
gen = sorted(set().union(*(carriers[k] for k, _ in CASS)),
             key=lambda a: (a != "MF6", org.get(a, ""), a))
M = np.array([[1 if g in carriers[k] else 0 for g in gen] for k, _ in CASS])
axB.imshow(M, cmap=plt.matplotlib.colors.ListedColormap([ABSENT, PRESENT]),
           aspect="auto", vmin=0, vmax=1, interpolation="nearest")


def short(g):
    """Species plus accession -- several carriers share a species name, so the
    species alone would label three different genomes identically."""
    s = (org.get(g, g).replace("Burkholderia ", "B. ")
         .replace("Paraburkholderia ", "P. "))
    s = s if len(s) <= 22 else s[:21] + "…"
    return f"{s}  {g}" if g != "MF6" else "B. sola MF6"


axB.set_xticks(range(len(gen)))
axB.set_xticklabels([short(g) for g in gen], rotation=90, fontsize=5.0)
for lbl, g in zip(axB.get_xticklabels(), gen):
    if g == "MF6":
        lbl.set_color(STAR); lbl.set_fontweight("bold")
axB.set_yticks(range(len(CASS)))
axB.set_yticklabels([lab for _, lab in CASS], fontsize=5.9)
axB.axhline(1.5, color="white", lw=2.2)   # 2 rows per cassette now, not 3
for sp in axB.spines.values():
    sp.set_visible(False)
axB.tick_params(length=0)
axB.set_title("B   Two RHS cassettes on pC3, different carrier sets\n"
              "      (filled = present at ≥60% id / ≥70% cov)",
              loc="left", fontweight="bold", fontsize=8.4)

# ---------------------------------------------------------------- panel C
axC = fig.add_subplot(gs[1, 1])
ys, labs, zs, cols, notes, hatches = [], [], [], [], [], []
# The two _Cterm loci were run by toxin12_phylo_test.py as a CONTRAST against
# their full-length query and were deliberately held out of the BH family (see
# its "[C-term contrast, not in BH]" tag), so they carry a raw p and no q. They
# must not be painted in the "not significant" grey -- they were never tested at
# that level. They get the same grey with a hatch and a raw-p label, and their
# own legend entry, so an uncorrected bar can never be read as a negative result.
for i, locus in enumerate(ORDER):
    p = phylo.get(locus, {})
    z = p.get("z_pc3", "")
    labs.append(ring_label(locus))
    ys.append(i)
    if z == "":
        zs.append(np.nan); cols.append(NOTSIG); notes.append("n/a")
        hatches.append(None)
    elif not p.get("p_pc3_bh"):
        zs.append(float(z)); cols.append(NOTSIG)
        notes.append("p=" + p["p_pc3"] + " raw")
        hatches.append("///")
    else:
        zs.append(float(z))
        sig = p.get("significant_bh_0.05") == "True"
        cols.append(CLUSTER if sig else NOTSIG)
        notes.append("q=" + p["p_pc3_bh"])
        hatches.append(None)
axC.axvline(0, color="#888888", lw=0.7)
# Two loci are strongly OVER-dispersed (z = +32, +39) because their carriers are
# mostly pC3-negative. Plotting them to scale compresses the informative range
# (-5..0) to nothing, so the axis is clipped and they are drawn as arrows with
# their true value written on.
XLIM = (-7.6, 6.0)
for y, z, c, n, hh in zip(ys, zs, cols, notes, hatches):
    if not np.isfinite(z):
        axC.text(-7.3, y, f"no pC3 null: more carriers than the\n"
                          f"{_n_pc3_tips} pC3-positive tips",
                 fontsize=4.9, va="center", color="#777777")
        continue
    zc = max(min(z, XLIM[1]), XLIM[0])
    axC.barh([y], [zc], color=c, height=0.6, hatch=hh,
             edgecolor="#5A6068" if hh else "none", linewidth=0.4)
    if z > XLIM[1]:
        axC.annotate("", xy=(XLIM[1] * 0.99, y), xytext=(XLIM[1] * 0.80, y),
                     arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8))
        axC.text(XLIM[1] * 0.76, y, f"z = +{z:.0f}", fontsize=5.0,
                 va="center", ha="right", color="#333333")
    else:
        axC.text(z - 0.35, y, n, fontsize=5.0, va="center", ha="right",
                 color="#333333")
axC.set_xlim(*XLIM)
axC.set_yticks(ys); axC.set_yticklabels(labs, fontsize=5.7)
# extra blank rows at the bottom: the legend gained a third entry and at the old
# limit it sat on top of the "no pC3 null" annotations
axC.set_ylim(len(ys) + 1.5, -1.5)
axC.set_xlabel("z of mean patristic distance\n(null = pC3-positive tips)",
               fontsize=6.4)
axC.tick_params(axis="x", labelsize=6)
axC.set_title("C   Clustered on the tree?\n      ← tighter than chance   ·   "
              "looser than chance →", loc="left", fontweight="bold",
              fontsize=8.4)
axC.legend(handles=[Patch(facecolor=CLUSTER, label="clustered (BH q < 0.05)"),
                    Patch(facecolor=NOTSIG, label="not significant"),
                    Patch(facecolor=NOTSIG, hatch="///", edgecolor="#5A6068",
                          linewidth=0.4,
                          label="C-term contrast: raw p, not in the BH family")],
           frameon=False, fontsize=5.6, loc="lower right",
           bbox_to_anchor=(1.0, 0.0))

fig.suptitle("MF6 carries two independent RHS toxin cassettes on pC3: one "
             "confined to a single clade,\nthe other scattered across the "
             "Burkholderia cepacia complex", fontsize=10, fontweight="bold",
             y=0.988)

# isolation-source key: a full-width strip at the bottom, as in Figure 1. It
# cannot go in axK -- twelve more swatches there collided with panel C. Only the
# categories actually present on these tips are listed, in the fixed
# vocabulary order so the two figures' legends read the same way.
_present_hosts = [h for h in figstyle.HOST_COLORS
                  if any(host.get(t.name, "unknown") == h for t in tips)]
_leg = fig.legend(handles=[Patch(facecolor=figstyle.HOST_COLORS[h],
                                 edgecolor="#999999", linewidth=0.3,
                                 label=h.replace("_", " "))
                           for h in _present_hosts],
                  loc="lower center", ncol=6, frameon=False, fontsize=6.3,
                  handlelength=1.2, columnspacing=1.1, labelspacing=0.35,
                  bbox_to_anchor=(0.5, 0.008),
                  title="isolation source (ring 2) — same categories as Figure 1",
                  title_fontsize=7)
_leg.get_title().set_fontweight("bold")
_leg.set_zorder(10)

ok = figstyle.save(fig, str(FIG / "fig11_toxin_tree"))
print("OK" if ok else "FAILED vector-text check")
