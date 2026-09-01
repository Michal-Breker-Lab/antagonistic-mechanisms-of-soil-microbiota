"""Minimal radial tree layout for annotated circular phylograms.

Written rather than pulled from ete3/toytree because neither is installed on the
Drive machines and both are heavy dependencies for what is ~80 lines of geometry.
Returns plain coordinates so the caller controls every drawing decision.
"""
import numpy as np


class RadialLayout:
    """Radial layout, optionally with whole clades collapsed to a single wedge.

    `collapse` is an iterable of Bio.Phylo clades to render as one triangular
    wedge instead of one slot per descendant tip. Each collapsed clade consumes
    `collapse_weight` tip-equivalents of arc rather than len(clade), which is the
    point: an outgroup of 125 tips can be given the arc of 12 and hand the rest
    back to the ingroup. Descendants of a collapsed clade get no coordinates, so
    nothing inside it is drawn; its depth is still measured from the real branch
    lengths, so the wedge is not a lie about how divergent it is.

    With `collapse=()` the layout is the previous one rotated by half a tip
    pitch: slots now have real angular width, so a tip sits at the centre of its
    slot rather than at its leading edge. That is what lets the rings use each
    slot's own half-width instead of one global step.
    """

    def __init__(self, tree, start_deg=90.0, extent_deg=350.0, use_branch_lengths=True,
                 collapse=(), collapse_weight=12.0):
        """tree: a Bio.Phylo tree (will be ladderized by the caller if wanted).

        `collapse` may be a sequence of clades -- all given `collapse_weight`
        tip-equivalents of arc -- or a {clade: weight} mapping when different
        wedges deserve different amounts of room (a 125-tip outgroup needs more
        lettering space than a 25-tip near-clone block).
        """
        self.tree = tree
        self.all_tips = tree.get_terminals()
        if hasattr(collapse, "items"):
            self.collapsed_clades = list(collapse.keys())
            weight_of = {id(c): float(w) for c, w in collapse.items()}
        else:
            self.collapsed_clades = list(collapse)
            weight_of = {id(c): float(collapse_weight) for c in self.collapsed_clades}
        collapse_ids = set(weight_of)

        # --- slot order: walk in tip order, but stop at a collapsed clade -----
        slots = []            # (kind, obj, weight); kind in {"tip", "clade"}
        def _walk(cl):
            if id(cl) in collapse_ids:
                slots.append(("clade", cl, weight_of[id(cl)]))
                return
            if not cl.clades:
                slots.append(("tip", cl, 1.0))
                return
            for c in cl.clades:
                _walk(c)
        _walk(tree.root)

        self.tips = [o for k, o, _ in slots if k == "tip"]
        n_slots = len(slots)
        span = np.deg2rad(extent_deg)
        off = np.deg2rad(start_deg)
        total_w = sum(w for _, _, w in slots) or 1.0
        # One slot would otherwise sit at angle 0 with no width; keep the old
        # convention that slot centres span [0, extent] inclusive.
        self.angle = {}
        self.halfwidth = {}
        self.wedge_span = {}
        cum = 0.0
        for kind, obj, w in slots:
            a0 = off + span * cum / total_w
            a1 = off + span * (cum + w) / total_w
            self.angle[id(obj)] = 0.5 * (a0 + a1)
            self.halfwidth[id(obj)] = 0.5 * (a1 - a0)
            if kind == "clade":
                self.wedge_span[id(obj)] = (a0, a1)
            cum += w
        if n_slots == 1:
            only = slots[0][1]
            self.angle[id(only)] = off

        # internal node angle = mean of children (collapsed clades act as leaves)
        for cl in tree.get_nonterminals(order="postorder"):
            if id(cl) in collapse_ids or id(cl) in self.angle:
                continue
            a = [self.angle[id(c)] for c in cl.clades if id(c) in self.angle]
            self.angle[id(cl)] = float(np.mean(a)) if a else 0.0

        # radius = root-to-node distance; do not descend into collapsed clades
        self.radius = {}
        root = tree.root
        self.radius[id(root)] = 0.0
        for cl in tree.get_nonterminals(order="preorder"):
            if id(cl) not in self.radius or id(cl) in collapse_ids:
                # not reachable, or IS a collapsed clade -- either way do not
                # descend. Descending into a collapsed clade would give its
                # children radii but no angles, which is the shape of a crash.
                continue
            for c in cl.clades:
                bl = c.branch_length if (use_branch_lengths and c.branch_length) else 0.0
                if not use_branch_lengths:
                    bl = 1.0
                self.radius[id(c)] = self.radius[id(cl)] + max(bl, 0.0)
        # true crown depth of each collapsed clade, in the same pre-normalised units
        self._crown = {}
        for cl in self.collapsed_clades:
            base = self.radius[id(cl)]
            deepest = 0.0
            for tp in cl.get_terminals():
                deepest = max(deepest, cl.distance(tp))
            self._crown[id(cl)] = base + deepest

        # Normalise on the DRAWN tips only. If a deep outgroup were allowed into
        # this maximum it would set the scale and re-compress the ingroup -- the
        # exact problem collapsing it is meant to solve. Collapsed crowns are
        # normalised by the same divisor, so a wedge deeper than the ingroup
        # simply reports a value > 1 and the caller caps it with a break marker
        # rather than silently rescaling everything else.
        drawn = [self.radius[id(t)] for t in self.tips] or [1.0]
        self.max_r = max(drawn) or 1.0
        for k in self.radius:
            self.radius[k] /= self.max_r
        for k in self._crown:
            self._crown[k] /= self.max_r

    def collapsed_wedges(self):
        """-> [(clade, a0, a1, r_stem, r_crown)] for the caller to draw."""
        out = []
        for cl in self.collapsed_clades:
            a0, a1 = self.wedge_span[id(cl)]
            out.append((cl, a0, a1, self.radius[id(cl)], self._crown[id(cl)]))
        return out

    def segments(self):
        """-> (radial_lines, arcs) ready for LineCollection."""
        radial, arcs = [], []
        for cl in self.tree.get_nonterminals(order="preorder"):
            if id(cl) not in self.radius:
                continue          # inside a collapsed clade -- nothing is drawn
            ra, aa = self.radius[id(cl)], self.angle[id(cl)]
            childs = [c for c in cl.clades if id(c) in self.radius]
            if not childs:
                continue
            angs = [self.angle[id(c)] for c in childs]
            arcs.append((ra, min(angs), max(angs)))
            for c in childs:
                rc, ac = self.radius[id(c)], self.angle[id(c)]
                radial.append(((ra, ac), (rc, ac)))
            del aa
        return radial, arcs

    def tip_angle(self, tip):
        return self.angle[id(tip)]

    def tip_radius(self, tip):
        return self.radius[id(tip)]


def polar_to_xy(r, a):
    return r * np.cos(a), r * np.sin(a)


def arc_points(r, a0, a1, n=40):
    a = np.linspace(a0, a1, max(n, 2))
    return r * np.cos(a), r * np.sin(a)


def draw_tree(ax, lay, lw=0.35, color="#333333"):
    from matplotlib.collections import LineCollection
    radial, arcs = lay.segments()
    segs = []
    for (r0, a0), (r1, a1) in radial:
        segs.append([polar_to_xy(r0, a0), polar_to_xy(r1, a1)])
    for r, a0, a1 in arcs:
        xs, ys = arc_points(r, a0, a1, n=max(4, int(np.rad2deg(a1 - a0))))
        segs.append(list(zip(xs, ys)))
    ax.add_collection(LineCollection(segs, colors=color, linewidths=lw, zorder=2))


def draw_ring(ax, lay, values, colors, r0, r1, default="#EEEEEE"):
    """values: dict tip_name -> key into colors. Draws one wedge per tip."""
    from matplotlib.patches import Wedge
    tips = lay.tips
    n = len(tips)
    if n < 2:
        return
    for t in tips:
        a = lay.tip_angle(t)
        hw = lay.halfwidth[id(t)]
        c = colors.get(values.get(t.name), default)
        ax.add_patch(Wedge((0, 0), r1, np.rad2deg(a - hw),
                           np.rad2deg(a + hw), width=r1 - r0,
                           facecolor=c, edgecolor="none", zorder=3))


def draw_bars(ax, lay, values, r0, maxlen, color="#1B3A6B", vmax=None):
    """radial bars scaled to maxlen (in axis units) beyond r0."""
    from matplotlib.patches import Wedge
    tips = lay.tips
    vals = [values.get(t.name, 0) or 0 for t in tips]
    vmax = vmax or (max(vals) if vals else 1) or 1
    for t in tips:
        v = values.get(t.name, 0) or 0
        if v <= 0:
            continue
        w = maxlen * v / vmax
        a = lay.tip_angle(t)
        hw = lay.halfwidth[id(t)] * 0.85
        ax.add_patch(Wedge((0, 0), r0 + w, np.rad2deg(a - hw),
                           np.rad2deg(a + hw), width=w,
                           facecolor=color, edgecolor="none", zorder=3))
    return vmax
