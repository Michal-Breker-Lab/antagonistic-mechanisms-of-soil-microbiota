"""Shared display conventions for the chromosome-1 tree figures.

Figures 4, 10 and 11 all draw the same 763-tip chromosome-1 tree. They used to
each midpoint-root it independently, which was both unreadable and inconsistent
with the rooting the report describes. This module holds the two decisions so
they cannot drift apart between figures.

Nothing here changes a topology or a branch length. Both operations are display
conventions: the rooting moves where the drawing starts, and the collapsing
hides tips that are already indistinguishable at the scale of the figure. No
statistic in the report is computed from these.

--- 1. Rooting -------------------------------------------------------------

Midpoint rooting lands on the branch subtending the 14 *Mycetohabitans*, so the
drawing was split 14 tips against 749; because radius is root-to-node distance,
the deep between-genera divergences then set the scale and the clade holding
every pC3 carrier got 32% of the radius for 84% of the tips.

`ingroup_root` instead roots on the ingroup stem. The ingroup is defined from the
TOPOLOGY, not from the organism labels, because the labels are polyphyletic
here: eight "Burkholderia sp." genomes sit inside the
Paraburkholderia/Caballeronia assemblage, so "all non-Burkholderia tips" is not
a clade and cannot be rooted on. The definition used is the largest clade that
contains every pC3 carrier and no non-Burkholderia-labelled tip. On this tree it
is unique and unambiguous -- 638 tips, 100% UFBoot -- and its parent is the node
that first admits Paraburkholderia.

--- 2. Near-clone collapsing ----------------------------------------------

Without dereplication (D11) the tree is mostly near-identical genomes: the
B. pseudomallei/mallei block alone is 232 tips whose entire crown spans 0.0055
subs/site. Drawn one tip per slot they produce a dense uninformative rim and a
hollow interior, and they crowd out the structure that carries the result.

`nearclone_clades` selects clades to collapse on three conditions:

    * crown height < MAX_CROWN  -- genuinely near-identical, not merely large
    * >= MIN_SIZE tips          -- collapsing 3 tips buys nothing
    * every tip pC3-NEGATIVE    -- never hide a carrier

The third condition is the important one and is not arbitrary: pC3 presence is
the subject of these figures, so a carrier must always be individually visible.
It also means a wedge's pC3 ring cell is a real measured value rather than a
"mixed" placeholder. A consequence worth stating in captions: collapsing can
only ever hide pC3-negative genomes, so it cannot flatter the pC3 signal.
"""
from __future__ import annotations

import collections
import csv
from pathlib import Path

import numpy as np

MAX_CROWN = 0.008      # subs/site; B. pseudomallei/mallei crown is 0.0055
MIN_SIZE = 20          # tips


def genus_of(acc, org):
    o = org.get(acc, "") or ""
    return o.split()[0] if o else "?"


def crown_height(clade):
    """Deepest root-to-tip distance inside `clade`, in substitutions/site."""
    return max(clade.distance(t) for t in clade.get_terminals())


def ingroup_root(tree, c3_present, org, verbose=True):
    """Root `tree` in place on the ingroup stem.

    Returns (ingroup_clade, outgroup_clade). `tree` must already be rooted
    somehow (Bio.Phylo needs a root to walk from); midpoint is fine.
    """
    carriers = {t.name for t in tree.get_terminals() if c3_present.get(t.name)}
    if not carriers:
        raise ValueError("no pC3 carriers on this tree - cannot define an ingroup")
    cands = [cl for cl in tree.get_nonterminals()
             if carriers <= {t.name for t in cl.get_terminals()}
             and all(genus_of(t.name, org) == "Burkholderia"
                     for t in cl.get_terminals())]
    if not cands:
        raise ValueError("no all-Burkholderia clade contains every pC3 carrier")
    ingroup = max(cands, key=lambda c: len(c.get_terminals()))
    names = sorted(t.name for t in ingroup.get_terminals())

    # Bio.Phylo roots AT the MRCA of the targets, so the result is a
    # trifurcation: outgroup plus the ingroup's two daughter lineages. That is
    # deliberate -- it spends no radius on an ingroup stem nobody reads.
    tree.root_with_outgroup(*[{"name": n} for n in names])
    tree.ladderize()

    inset = set(names)
    og = [c for c in tree.root.clades
          if not ({t.name for t in c.get_terminals()} & inset)]
    if len(og) != 1:
        raise ValueError(f"expected one outgroup clade off the root, got {len(og)}")
    outgroup = og[0]
    if verbose:
        comp = collections.Counter(genus_of(t.name, org)
                                   for t in outgroup.get_terminals())
        print(f"  ingroup: {len(names)} tips (all Burkholderia, all "
              f"{len(carriers)} carriers)")
        print(f"  outgroup: {len(outgroup.get_terminals())} tips -> "
              + ", ".join(f"{k} {v}" for k, v in comp.most_common()))
        n_pos = sum(1 for t in outgroup.get_terminals() if c3_present.get(t.name))
        assert n_pos == 0, "a pC3 carrier ended up in the outgroup"
    return ingroup, outgroup


def nearclone_clades(tree, c3_present, org=None, skip=(), max_crown=MAX_CROWN,
                     min_size=MIN_SIZE, verbose=True):
    """Maximal clades to collapse. `skip` = clades to stay out of (the outgroup)."""
    skip_ids = set()
    for cl in skip:
        skip_ids.update(id(x) for x in cl.find_clades())
    kept, done = [], set()
    for cl in tree.find_clades(order="preorder"):
        if cl.is_terminal() or id(cl) in done or id(cl) in skip_ids:
            continue
        tl = cl.get_terminals()
        if len(tl) < min_size:
            continue
        if any(c3_present.get(t.name) for t in tl):
            continue                       # never collapse a pC3 carrier
        if crown_height(cl) >= max_crown:
            continue
        kept.append(cl)
        done.update(id(x) for x in cl.find_clades())
    if verbose:
        hidden = sum(len(c.get_terminals()) for c in kept)
        print(f"  near-clone wedges: {len(kept)} hiding {hidden} tips "
              f"(crown < {max_crown}, n >= {min_size}, all pC3-negative)")
        for c in sorted(kept, key=lambda x: -len(x.get_terminals())):
            print(f"     {len(c.get_terminals()):>4} tips  "
                  f"crown={crown_height(c):.5f}  {clade_label(c, org or {})}")
    return kept


def clade_label(clade, org, italic=False, abbrev=False):
    """Dominant 'Genus species' label, with the clade size.

    `abbrev` shortens the genus to an initial, which is what makes a wedge label
    fit beside a 763-tip circular tree. The size given is always the FULL clade
    size, never the dominant species' count -- a wedge labelled (232) contains
    232 genomes even though 35 of them carry a different species name (here,
    B. mallei inside B. pseudomallei, which is expected: mallei is a clonal
    derivative of pseudomallei, not a separate lineage).
    """
    names = [" ".join((org.get(t.name, "") or "?").split()[:2])
             for t in clade.get_terminals()]
    lab, n = collections.Counter(names).most_common(1)[0]
    tot = len(names)
    if abbrev and " " in lab:
        g, sp = lab.split(None, 1)
        lab = f"{g[0]}. {sp}"
    if italic and lab != "?":
        lab = "$\\it{" + "\\ ".join(lab.split()) + "}$"
    return f"{lab} ({tot})"


def draw_wedge(ax, rt, a0, a1, r_stem, r_cap, facecolor="#E8E4DC",
               edgecolor="#7A7A7A", lw=0.6, zorder=2, break_marker=True):
    """Filled triangle for a collapsed clade, optionally with a scale break.

    The break is drawn whenever the clade's true depth exceeds what is shown, so
    a truncated wedge never silently understates divergence.
    """
    from matplotlib.patches import Polygon
    xs, ys = rt.arc_points(r_cap, a0, a1, n=24)
    pts = [rt.polar_to_xy(r_stem, 0.5 * (a0 + a1))] + list(zip(xs, ys))
    ax.add_patch(Polygon(pts, closed=True, facecolor=facecolor,
                         edgecolor=edgecolor, linewidth=lw, zorder=zorder))
    if break_marker:
        for frac in (0.80, 0.835):
            r = r_stem + (r_cap - r_stem) * frac
            bx, by = rt.arc_points(r, a0, a1, n=12)
            ax.plot(bx, by, color="#FFFFFF", lw=1.5, zorder=zorder + 1,
                    solid_capstyle="butt")
            ax.plot(bx, by, color=edgecolor, lw=0.5, zorder=zorder + 1,
                    solid_capstyle="butt")


def ring_cell(ax, a0, a1, r0, r1, facecolor, zorder=3):
    """One ring cell spanning a collapsed wedge's arc."""
    from matplotlib.patches import Wedge
    ax.add_patch(Wedge((0, 0), r1, np.rad2deg(a0), np.rad2deg(a1),
                       width=r1 - r0, facecolor=facecolor, edgecolor="none",
                       zorder=zorder))


# --- 3. Isogenic lab derivatives ------------------------------------------
# The 791-genome bac120 set was built from "every complete/chromosome-level
# assembly with isolation metadata", and that filter cannot tell an isolate from
# a laboratory derivative of one. It let in BioProject PRJEB40633 (LM-UGent):
# one parent strain, B. sola R-12632 from a maize rhizosphere in Italy, 1996,
# plus 37 TRANSPOSON MUTANTS of that same strain, sequenced to map the genes
# behind its antibacterial activity (Depoorter, Coenye & Vandamme 2021, Appl
# Environ Microbiol 87:e01169-21). Every mutant inherits the parent BioSample's
# metadata verbatim, so they present as 38 independent maize-rhizosphere
# isolates. They are one genome.
#
# Left in, they are not merely redundant, they are misleading exactly where MF6
# sits: 34 of the 38 share one identical bac120 sequence, all 38 are pC3+ and 37
# carry the full ten-locus toxin set, so any "N of M tips carry X" count in
# MF6's neighbourhood is inflated ~8-fold by a single strain.
#
# A whole-set scan of NCBI BioProject/BioSample records (strain, submitter id,
# sample name, project title matched against mutant/derivative/transposon/
# knockout/substr.) found NO other lab-derivative block in the 790 assemblies --
# this is the only one.
#
# The parent is KEPT. This is a display/redundancy decision, not a re-inference:
# the tree and the alignment are untouched, the mutants are simply not drawn.
DERIVATIVES_TSV = (Path(__file__).resolve().parents[2]
                   / "bac120" / "tables" / "lab_derivatives_excluded.tsv")


def lab_derivatives(path=DERIVATIVES_TSV):
    """Accessions that are isogenic laboratory derivatives of another genome."""
    if not path.exists():
        raise SystemExit(f"missing {path} - run the derivative scan first")
    with open(path, newline="") as fh:
        return {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}


def drop_lab_derivatives(tree, path=DERIVATIVES_TSV, verbose=True):
    """Prune isogenic lab derivatives; return the accessions actually removed.

    Pruning happens BEFORE any layout or count, so no figure statistic is
    computed over them. Parents are asserted to survive: dropping a mutant block
    while also losing its parent would delete a real genome from the panel.
    """
    tab = lab_derivatives(path)
    on_tree = {t.name for t in tree.get_terminals()}
    drop = [a for a in tab if a in on_tree]
    for acc in drop:
        tree.prune(target=acc)
    left = {t.name for t in tree.get_terminals()}
    parents = {r["parent_accession"] for a, r in tab.items() if a in drop}
    missing = sorted(p for p in parents if p not in left)
    assert not missing, f"parent strain(s) not on the pruned tree: {missing}"
    if verbose:
        by_parent = collections.Counter(tab[a]["parent_accession"] for a in drop)
        for p, n in by_parent.most_common():
            print(f"  lab derivatives dropped: {n} isogenic mutants of {p} "
                  f"({tab[drop[0]]['parent_strain']}) - kept the parent")
    return drop
