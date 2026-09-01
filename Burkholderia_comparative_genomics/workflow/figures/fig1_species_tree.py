#!/usr/bin/env python3
"""Figure 4 - chromosome-1 species tree with c3 presence, host and compartment.

The tree is built from chromosome-1 core genes ONLY, so c3 presence cannot have
influenced the topology it is being mapped onto.

Branches are painted by host. The painting rule is deliberately model-free: a
branch takes a host colour only if EVERY tip descended from it shares that host,
otherwise it is drawn grey. No ancestral-state reconstruction is performed and
none is implied -- grey means "descendants disagree", not "unknown ancestor".
Using parsimony or ML reconstruction here would assert ancestral habitats that
these data cannot support.

Rings, inner to outer:
  1. c3 presence / absence
  2. host category

The plant-compartment ring and the c3-size bars were removed 2026-08-09 at
Moshe's request -- the compartment ring was a sub-annotation of a single host
category (so mostly empty), and c3 size is already carried by figure 5.
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import figstyle  # noqa: E402
import radialtree as rt  # noqa: E402
import tree_display as td  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from Bio import Phylo  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

D = Path(__file__).resolve().parent.parent
TAB, FIG = D / "tables", D / "figures"
FIG.mkdir(exist_ok=True)
TREE = TAB / "chr1_core.treefile"

if not TREE.exists():
    sys.exit(f"tree not found: {TREE} - run Stage 8 first")

tree = Phylo.read(str(TREE), "newick")
# IQ-TREE writes an unrooted tree, so a rooting has to be chosen for display.
#
# This figure used to be MIDPOINT-rooted, and that was the wrong choice for
# reading it. Midpoint rooting lands on the branch subtending the 14
# Mycetohabitans, so the drawing was split 14 tips against 749, and because the
# radius is root-to-node distance the deep between-genera divergences ate the
# scale: the *Burkholderia* sensu stricto clade -- where every pC3 carrier is --
# got 32% of the radius to hold 84% of the tips.
#
# Instead the tree is rooted on the INGROUP STEM and the outgroup is collapsed to
# one labelled wedge (see below). Rooting is still a display convention and no
# statistic in the report depends on it, but this one is also better supported
# than midpoint: the ingroup/outgroup split carries 100% UFBoot.
tree.root_at_midpoint()
tree.ladderize()
tips_all = tree.get_terminals()
print(f"tips: {len(tips_all)}")

# ---------------- trait tables ----------------
host, org = {}, {}
for r in csv.DictReader(open(TAB / "host_categories.tsv"), delimiter="\t"):
    host[r["accession"]] = r["host_category"]
    org[r["accession"]] = r["organism_name"]

c3_present = {}
for r in csv.DictReader(open(TAB / "secondary_replicon_clusters.tsv"), delimiter="\t"):
    if r["is_c3"] == "True":
        c3_present[r["accession"]] = True
for t in tips_all:
    c3_present.setdefault(t.name, False)

# ---------------- ingroup / outgroup / near-clone wedges ----------------
# Both display conventions live in tree_display so figures 4, 10 and 11 cannot
# drift apart. See that module for why the ingroup is defined from the topology
# rather than the (polyphyletic) organism labels, and why a collapsed clade is
# never allowed to contain a pC3 carrier.
print("display conventions:")
ingroup, outgroup = td.ingroup_root(tree, c3_present, org)
out_tips = outgroup.get_terminals()
nearclones = td.nearclone_clades(tree, c3_present, org, skip=[outgroup])
_hidden = {t.name for cl in nearclones for t in cl.get_terminals()}
tips = [t for t in tree.get_terminals()
        if t.name not in {x.name for x in out_tips} and t.name not in _hidden]

# MF6 is a lab isolate with no NCBI BioSample, so it is carried as a manually
# curated row in host_categories.tsv rather than hard-coded here -- that file is
# the single source of truth shared with fig6, which previously disagreed and
# drew MF6 as "unknown".
if "MF6" not in host:
    print("WARNING: MF6 missing from host_categories.tsv - it will draw as unknown")

n_pres = sum(1 for t in tips_all if c3_present.get(t.name))
print(f"c3 present on tree: {n_pres}/{len(tips_all)}")

print(f"MF6 host: {host.get('MF6', 'MISSING')}  organism: {org.get('MF6', 'MISSING')}")

# ---------------- draw ----------------
fig = plt.figure(figsize=(7.8, 8.2))
ax = fig.add_axes([0.02, 0.095, 0.96, 0.825])
ax.set_aspect("equal"); ax.axis("off")

# Weights are in TIP-EQUIVALENTS of arc. The outgroup gets more than the
# near-clone blocks because it carries three lines of lettering.
OG_WEIGHT, NC_WEIGHT = 22.0, 10.0
_collapse = {outgroup: OG_WEIGHT}
_collapse.update({cl: NC_WEIGHT for cl in nearclones})
lay = rt.RadialLayout(tree, start_deg=92, extent_deg=352, collapse=_collapse)
_slots = len(tips) + OG_WEIGHT + NC_WEIGHT * len(nearclones)
print(f"drawn elements: {len(tips)} tips + {len(nearclones)} near-clone wedges "
      f"+ 1 outgroup wedge (from {len(tips_all)} tips)")
print(f"angular pitch: {352.0/_slots:.4f} deg/slot "
      f"(was {352.0/len(tips_all):.4f} over all {len(tips_all)} tips, "
      f"{_slots and 352.0/_slots/(352.0/len(tips_all)):.2f}x)")

# ---- branch painting by host ----------------------------------------------
# A node's state is its host if all descendant tips agree, else None (grey).
state = {}
for t in tree.get_terminals():
    state[id(t)] = host.get(t.name, "unknown")
for cl in tree.get_nonterminals(order="postorder"):
    s = {state[id(c)] for c in cl.clades}
    state[id(cl)] = s.pop() if len(s) == 1 else None

GREY = "#C8C8C8"


def colour_of(node):
    st = state.get(id(node))
    return figstyle.HOST_COLORS.get(st, GREY) if st else GREY


segs, cols = [], []
for cl in tree.get_nonterminals(order="preorder"):
    if id(cl) not in lay.radius:
        continue                      # inside the collapsed outgroup
    ra = lay.radius[id(cl)]
    # the connecting arc sits at the parent's radius -> parent's state
    angs = [lay.angle[id(c)] for c in cl.clades if id(c) in lay.radius]
    if len(angs) > 1:
        xs, ys = rt.arc_points(ra, min(angs), max(angs),
                               n=max(4, int(np.rad2deg(max(angs) - min(angs)))))
        segs.append(list(zip(xs, ys))); cols.append(colour_of(cl))
    # each radial branch is coloured by the CHILD it leads to
    for c in cl.clades:
        if id(c) not in lay.radius:
            continue
        rc, ac = lay.radius[id(c)], lay.angle[id(c)]
        segs.append([rt.polar_to_xy(ra, ac), rt.polar_to_xy(rc, ac)])
        cols.append(colour_of(c))
ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.55, zorder=2))
n_painted = sum(1 for c in cols if c != GREY)
print(f"branch segments painted by host: {n_painted}/{len(cols)} "
      f"({100*n_painted/len(cols):.0f}%)")

R0 = 1.04
# dotted extensions from each tip out to the rings -- in a phylogram tips end at
# different radii, and without these the rings float free of the tree
_ext = [[rt.polar_to_xy(lay.tip_radius(t), lay.tip_angle(t)),
         rt.polar_to_xy(R0, lay.tip_angle(t))] for t in tips]
ax.add_collection(LineCollection(_ext, colors="#BBBBBB", linewidths=0.22,
                                 linestyles=(0, (1, 2)), zorder=1))

rt.draw_ring(ax, lay, {t.name: c3_present.get(t.name, False) for t in tips},
             figstyle.C3_COLORS, R0, R0 + 0.055, default="#D9D9D9")
rt.draw_ring(ax, lay, {t.name: host.get(t.name, "unknown") for t in tips},
             figstyle.HOST_COLORS, R0 + 0.068, R0 + 0.123, default="#EEEEEE")

# ---- the collapsed wedges --------------------------------------------------
# A wedge is drawn at its TRUE depth whenever that fits inside the ingroup
# radius, which is the case for every near-clone block. Only the outgroup
# exceeds it (2.6x), and that one alone gets a break marker plus the real number
# in its label -- capping it silently would understate how divergent those
# genera are, and letting it set the scale is what made this figure unreadable.
OG_CAP = 0.86
for cl, a0, a1, r_stem, r_crown in lay.collapsed_wedges():
    truncated = r_crown > 1.0
    r_cap = OG_CAP if truncated else r_crown
    td.draw_wedge(ax, rt, a0, a1, r_stem, r_cap, break_marker=truncated)
    # dotted leaders across the wedge, exactly as for individual tips: a
    # near-clone block is a radial SLIVER (232 genomes spanning 0.0055
    # subs/site), so without these the eye cannot connect it to its ring cells.
    _la = np.linspace(a0, a1, 7)[1:-1]
    ax.add_collection(LineCollection(
        [[rt.polar_to_xy(r_cap, _a), rt.polar_to_xy(R0, _a)] for _a in _la],
        colors="#BBBBBB", linewidths=0.22, linestyles=(0, (1, 2)), zorder=1))

    # Rings continue across every wedge. The pC3 cell is a measured value, not a
    # placeholder: nothing with a carrier in it is ever collapsed.
    td.ring_cell(ax, a0, a1, R0, R0 + 0.055, figstyle.C3_COLORS[False])
    _h = {host.get(t.name, "unknown") for t in cl.get_terminals()}
    # Mixed hosts inside a wedge reuse GREY -- the same colour the branch
    # painting uses for "descendants disagree". Same meaning, same colour.
    _hc = figstyle.HOST_COLORS.get(next(iter(_h)), "#EEEEEE") if len(_h) == 1 \
        else GREY
    td.ring_cell(ax, a0, a1, R0 + 0.068, R0 + 0.123, _hc)

    am = 0.5 * (a0 + a1)
    if cl is outgroup:
        lx, ly = rt.polar_to_xy(R0 + 0.215, am)
        ax.text(lx, ly, f"outgroup, {len(cl.get_terminals())} genomes\n"
                        f"$\\it{{Paraburkholderia}}$ 84 · $\\it{{Caballeronia}}$ 16 · "
                        f"$\\it{{Mycetohabitans}}$ 14 · $\\it{{Trinickia}}$ 3\n"
                        f"collapsed; true depth {r_crown:.1f}× this radius · no pC3",
                fontsize=5.9, color="#333333",
                # extend AWAY from the neighbouring near-clone labels rather
                # than centring on the wedge, which collided with them
                ha="right" if np.cos(am) < 0 else "left",
                va="bottom" if np.sin(am) > 0 else "top", zorder=9,
                linespacing=1.35)
    else:
        deg = np.rad2deg(am) % 360.0
        flip = 90.0 < deg < 270.0
        lx, ly = rt.polar_to_xy(R0 + 0.135, am)
        ax.text(lx, ly, td.clade_label(cl, org, italic=True, abbrev=True),
                fontsize=6.0, color="#333333", zorder=9,
                rotation=deg + 180 if flip else deg,
                rotation_mode="anchor",
                ha="right" if flip else "left", va="center")

# mark MF6 (restored: this block must stay AFTER the wedges so the star draws
# on top of them, and it is the only tip label on the figure)
for t in tips:
    if t.name == "MF6":
        a = lay.tip_angle(t)
        x, y = rt.polar_to_xy(R0 + 0.185, a)
        ax.scatter([x], [y], s=64, marker="*", c="#D55E00",
                   edgecolors="black", linewidths=0.5, zorder=10)
        xt, yt = rt.polar_to_xy(R0 + 0.235, a)
        ax.text(xt, yt, "MF6", fontsize=7.5, fontweight="bold", color="#D55E00",
                ha="left" if np.cos(a) > 0 else "right", va="center", zorder=10)
        break
else:
    print("WARNING: MF6 not among the drawn tips - star not placed")

lim = R0 + 0.46
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

# scale bar (radii are normalised to root-to-deepest-tip distance)
depth = lay.max_r
_n = 10 ** np.floor(np.log10(depth / 4))
step = _n * max(1, round((depth / 4) / _n))
ax.plot([-lim + 0.04, -lim + 0.04 + step / depth], [-lim + 0.10, -lim + 0.10],
        lw=1.2, color="#000000")
ax.text(-lim + 0.04 + step / (2 * depth), -lim + 0.14, f"{step:g} subs/site",
        ha="center", va="bottom", fontsize=7)

fig.text(0.5, 0.985, "Replicon 3 (pC3) presence across the chromosome-1 "
         "phylogeny of $\\it{Burkholderia}$ sensu lato",
         ha="center", va="top", fontsize=10.5, fontweight="bold")
fig.text(0.5, 0.958, f"{len(tips_all)} complete genomes, not dereplicated (D11); "
         f"tree from chromosome-1 core genes only — pC3 present in {n_pres} "
         f"({100*n_pres/len(tips_all):.0f}%). Rooted on the ingroup stem (100% "
         f"UFBoot); the {len(out_tips)}-genome outgroup and {len(nearclones)} "
         f"near-clone blocks ({len(_hidden)} genomes) are collapsed to wedges.\n"
         f"Branches painted by host where all descendants agree, grey otherwise.",
         ha="center", va="top", fontsize=7.2, color="#444444")

h1 = [Patch(facecolor=figstyle.C3_COLORS[True], label="c3 present"),
      Patch(facecolor=figstyle.C3_COLORS[False], label="c3 absent")]
present_hosts = [h for h in figstyle.HOST_COLORS
                 if any(host.get(t.name) == h for t in tips_all)]
h2 = [Patch(facecolor=figstyle.HOST_COLORS[h], label=h.replace("_", " "))
      for h in present_hosts]
h4 = [Line2D([0], [0], color=GREY, lw=1.6, label="branch: hosts disagree"),
      Patch(facecolor="#E8E4DC", edgecolor="#7A7A7A",
            label="collapsed clade (pC3-negative)")]

leg1 = fig.legend(handles=h1 + h2 + h4, loc="lower center", ncol=7,
                  frameon=False, fontsize=6.5, handlelength=1.2,
                  columnspacing=1.1, labelspacing=0.32,
                  bbox_to_anchor=(0.5, 0.012),
                  title="c3 presence  ·  host / isolation source",
                  title_fontsize=7)
leg1.get_title().set_fontweight("bold")
leg1.set_zorder(10)
fig.add_artist(leg1)

print("\nverifying vector text:")
figstyle.save(fig, str(FIG / "fig1_species_tree_c3"))
print("done")
