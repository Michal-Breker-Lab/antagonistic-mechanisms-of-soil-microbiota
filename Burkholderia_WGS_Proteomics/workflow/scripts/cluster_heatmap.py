#!/usr/bin/env python3
"""Clustermap of the day-2 core proteome: KEGG enrichment, z-scored abundances and
COG enrichment in one figure.
"""

import sys
from types import SimpleNamespace
import os
import re
import textwrap

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "mathtext.default": "regular",
})
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
import seaborn as sns

KEGG_STOPS = ["#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#0d366b"]
COG_STOPS = ["#e9ddf4", "#cdb0e6", "#a878d1", "#8348b8", "#5f2d90", "#3d1a60"]
KEGG_CMAP = LinearSegmentedColormap.from_list("kegg_blue", KEGG_STOPS)
COG_CMAP = LinearSegmentedColormap.from_list("cog_purple", COG_STOPS)
INK = "#0b0b0b"
GRID = "#D6D6D6"
LIGHT_BOX = "#BEBEBE"

STRAIN_ORDER = ["MF6", "27D6", "34F7"]
BAR_H = 0.62


def parse_args():
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        indir=sm.params.indir,
        prefix=sm.params.prefix,
        outdir=None,
        fdr=float(sm.params.fdr),
        min_count=int(sm.params.min_count),
        max_terms=5,
        vlim=2.0,
        png=False,
    )
    return args


def wrap_label(s, width=28, maxlines=2):
    """Wrap a long term name onto up to `maxlines` lines (each kept under ~30
    chars) instead of truncating, without ever splitting a word (or hyphenated
    word); only overflow past the last allowed line is trimmed with an ellipsis."""
    s = str(s)
    lines = textwrap.wrap(s, width=width, break_long_words=False,
                          break_on_hyphens=False) or [s]
    if len(lines) > maxlines:
        lines = lines[:maxlines]
        lines[-1] = lines[-1].rstrip() + "…"
    return "\n".join(lines)


def strain_chlamy_info(sample_cols):
    """(strains, chlamy flags, multi-Chlamy?, strain palette, Chlamy palette),
    matching core_clustermap.py so the strip + legend colors are consistent."""
    strains = [c.split("_")[0] for c in sample_cols]
    chlamy = ["_Chlamy_" in c for c in sample_cols]
    multi_chlamy = len(set(chlamy)) > 1
    strain_pal = dict(zip(STRAIN_ORDER, sns.color_palette("Set2", 3)))
    chlamy_pal = {False: "#dddddd", True: "#1b7837"}
    return strains, chlamy, multi_chlamy, strain_pal, chlamy_pal


def load_clusters(indir, prefix):
    """Return (df, sample_cols, zscored matrix). Sample columns are everything
    to the right of the 'cluster' column; the heatmap z-scores per protein
    across them, reproducing the original clustermap's rows."""
    path = os.path.join(indir, f"{prefix}_clusters.tsv")
    df = pd.read_csv(path, sep="\t", comment="#")
    cols = list(df.columns)
    sample_cols = cols[cols.index("cluster") + 1:]
    mat = df.set_index("protein_id")[sample_cols].astype(float)
    z = mat.sub(mat.mean(axis=1), axis=0)
    sd = mat.std(axis=1).replace(0, np.nan)
    z = z.div(sd, axis=0).fillna(0.0)
    return df, sample_cols, z


def load_enrichment(indir, prefix, kind, fdr, min_count):
    path = os.path.join(indir, f"{prefix}_{kind}_enrichment.tsv")
    e = pd.read_csv(path, sep="\t", comment="#")
    sig = e[(e["FDR"] <= fdr) & (e["k"] >= min_count)].copy()
    sig["neglog10_fdr"] = -np.log10(sig["FDR"].clip(lower=1e-300))
    return sig.sort_values(["cluster", "FDR"])


def group_boundaries(sample_cols):
    """Column indices where the STRAIN/Chlamy group changes (drops the trailing
    _<rep> from each name), for the vertical separators + top color strip."""
    keys = [re.sub(r"_\d+$", "", c) for c in sample_cols]
    return [i for i in range(1, len(keys)) if keys[i] != keys[i - 1]], keys


def panel_scale(sig, max_terms):
    """Shared x/color scale for one enrichment panel across all clusters:
    (fold-enrichment max, -log10FDR Normalize, slot count = max terms in any
    cluster). Uses only the terms that will actually be drawn (top max_terms).

    Bar length is FOLD ENRICHMENT, not the gene ratio the upstream figure used.
    The gene ratio k/n scales with how common a term is, not with how enriched it
    is: the retired superset ko01130 drew the longest bar in cluster 4 at k/n =
    0.31 with a fold enrichment of 1.40 - barely above background - while
    C5-branched dibasic acid metabolism at fold 3.45 was almost invisible beside
    it.  Fold enrichment is also what Fisher's exact actually tests, and what the
    companion bubble figure puts on its x axis, so the two now agree."""
    if sig.empty:
        return 0.1, Normalize(0, 1), 1
    shown = sig.groupby("cluster", group_keys=False).head(max_terms)
    gmax = float(shown["fold_enrichment"].max())
    norm = Normalize(vmin=float(shown["neglog10_fdr"].min()),
                     vmax=float(shown["neglog10_fdr"].max()))
    slots = int(shown.groupby("cluster").size().max())
    return gmax, norm, max(slots, 1)


def draw_bar_panel(ax, terms, side, cmap, norm, gmax, slots, label_field,
                   show_xaxis):
    """Horizontal bars: length = fold enrichment, color = -log10(FDR). Bars grow
    OUTWARD, away from the heatmap: side='left' (KEGG) bases sit on the right
    (heatmap side) and grow left, labels outside on the left; side='right'
    (COG) bases sit on the left and grow right, labels outside on the right.
    Bars are vertically centered in the band."""
    ax.set_ylim(slots - 0.5, -0.5)
    xhi = 1.0 + (gmax - 1.0) * 1.04
    ax.set_xlim((xhi, 1.0) if side == "left" else (1.0, xhi))
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}\u00d7"))
    ax.grid(axis="x", ls=":", lw=0.8, color=GRID, zorder=0)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(LIGHT_BOX)
        s.set_linewidth(0.8)
    ax.tick_params(axis="y", length=0)
    if side == "right":
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("right")

    if show_xaxis:
        ax.tick_params(axis="x", length=2.5, labelsize=7, bottom=True,
                       labelbottom=True, top=False, labeltop=False)
        ax.set_xlabel("fold enrichment", fontsize=8, labelpad=2)
    else:
        ax.tick_params(axis="x", bottom=False, labelbottom=False,
                       top=False, labeltop=False)

    if terms.empty:
        ax.set_yticks([])
        return

    m = len(terms)
    start = (slots - m) / 2.0
    ys = [start + j for j in range(m)]
    ax.barh(ys, terms["fold_enrichment"].values - 1.0, left=1.0, height=BAR_H,
            color=[cmap(norm(v)) for v in terms["neglog10_fdr"]],
            edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_yticks(ys)
    ax.set_yticklabels([wrap_label(t) for t in terms[label_field]], fontsize=8)


def main():
    args = parse_args()
    outdir = args.outdir or args.indir
    os.makedirs(outdir, exist_ok=True)

    df, sample_cols, z = load_clusters(args.indir, args.prefix)
    kegg = load_enrichment(args.indir, args.prefix, "pathway", args.fdr, args.min_count)
    cog = load_enrichment(args.indir, args.prefix, "cog", args.fdr, args.min_count)
    cog = cog.assign(label=cog["category"].astype(str) + "  " + cog["description"].astype(str))

    n_clusters = int(df["cluster"].max())
    sizes = df["cluster"].value_counts().sort_index()
    bounds, _ = group_boundaries(sample_cols)

    kegg_gmax, kegg_norm, kegg_slots = panel_scale(kegg, args.max_terms)
    cog_gmax, cog_norm, cog_slots = panel_scale(cog, args.max_terms)

    cluster_pal = sns.color_palette("tab10", max(n_clusters, 3))[:n_clusters]

    fig = plt.figure(figsize=(16, 1.55 * n_clusters + 2.3))
    gs = GridSpec(n_clusters, 4, width_ratios=[1.0, 0.13, 2.13, 1.0],
                  wspace=0.035, hspace=0.08)
    fig.subplots_adjust(top=0.87, bottom=0.16, left=0.17, right=0.83)

    top_axes, bottom_heat = {}, None
    for i in range(n_clusters):
        cid = i + 1
        cl_genes = df.loc[df["cluster"] == cid, "protein_id"].tolist()

        ax_k = fig.add_subplot(gs[i, 0])
        draw_bar_panel(ax_k, kegg[kegg["cluster"] == cid].head(args.max_terms),
                       "left", KEGG_CMAP, kegg_norm, kegg_gmax, kegg_slots,
                       "description", show_xaxis=(i == n_clusters - 1))

        ax_c = fig.add_subplot(gs[i, 1])
        ax_c.set_xlim(0, 1); ax_c.set_ylim(0, 1)
        ax_c.set_facecolor(cluster_pal[cid - 1])
        ax_c.set_xticks([]); ax_c.set_yticks([])
        for s in ax_c.spines.values():
            s.set_color("black"); s.set_linewidth(1.0)
        ax_c.annotate(f"C{cid}", xy=(0.5, 0.5), xytext=(0, 6),
                      textcoords="offset points", ha="center", va="center",
                      fontsize=11, fontweight="bold", color="white",
                      path_effects=[pe.withStroke(linewidth=1.8, foreground="black")])
        ax_c.annotate(str(int(sizes[cid])), xy=(0.5, 0.5), xytext=(0, -7),
                      textcoords="offset points", ha="center", va="center",
                      fontsize=8, color="white",
                      path_effects=[pe.withStroke(linewidth=1.4, foreground="black")])

        ax_h = fig.add_subplot(gs[i, 2])
        sns.heatmap(z.loc[cl_genes], ax=ax_h, cmap="RdBu_r", center=0,
                    vmin=-args.vlim, vmax=args.vlim, cbar=False,
                    xticklabels=False, yticklabels=False)
        ax_h.set_ylabel(""); ax_h.set_xlabel("")
        for s in ax_h.spines.values():
            s.set_visible(True); s.set_color("black"); s.set_linewidth(1.0)
        for x in bounds:
            ax_h.axvline(x, color="black", lw=0.8)

        ax_g = fig.add_subplot(gs[i, 3])
        draw_bar_panel(ax_g, cog[cog["cluster"] == cid].head(args.max_terms),
                       "right", COG_CMAP, cog_norm, cog_gmax, cog_slots,
                       "label", show_xaxis=(i == n_clusters - 1))

        if i == 0:
            top_axes = {"kegg": ax_k, "heat": ax_h, "cog": ax_g}
        if i == n_clusters - 1:
            bottom_heat = ax_h

    fig.canvas.draw()
    add_headers(fig, top_axes, sample_cols, args, kegg_norm, cog_norm)
    add_bottom_strip(fig, bottom_heat, sample_cols, bounds)

    pdf = os.path.join(outdir, f"{args.prefix}_heatmap_enrichment.pdf")
    fig.savefig(pdf, bbox_inches="tight")
    print(f"wrote {pdf}")
    svg = os.path.join(outdir, f"{args.prefix}_heatmap_enrichment.svg")
    fig.savefig(svg, bbox_inches="tight")
    print(f"wrote {svg}")
    if args.png:
        png = os.path.join(outdir, f"{args.prefix}_heatmap_enrichment.png")
        fig.savefig(png, bbox_inches="tight", dpi=150)
        print(f"wrote {png}")
    plt.close(fig)


def add_headers(fig, top_axes, sample_cols, args, kegg_norm, cog_norm):
    """Top of the figure: three aligned colorbars (KEGG blue, Expression RdBu,
    COG purple) hugging their panel tops, the three panel titles, and a single
    horizontal strain/Chlamy legend centered above them. The sample names now
    live at the bottom (add_bottom_strip), so all three headers align."""
    hm = top_axes["heat"].get_position()
    kp = top_axes["kegg"].get_position()
    cp = top_axes["cog"].get_position()
    kegg_cx = (kp.x0 + kp.x1) / 2
    heat_cx = (hm.x0 + hm.x1) / 2
    cog_cx = (cp.x0 + cp.x1) / 2

    cb_w, cb_h = 0.12, 0.011
    cb_y = hm.y1 + 0.008
    title_y = cb_y + cb_h + 0.030

    def cbar(cx, norm, cmap, label, ticks=None):
        cax = fig.add_axes([cx - cb_w / 2, cb_y, cb_w, cb_h])
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                          cax=cax, orientation="horizontal", ticks=ticks)
        cb.ax.xaxis.set_ticks_position("top")
        cb.ax.xaxis.set_label_position("top")
        cb.ax.tick_params(labelsize=7, pad=1)
        cb.set_label(label, fontsize=8, labelpad=3)
        if getattr(cb, "solids", None) is not None:
            cb.solids.set_rasterized(False)

    cbar(kegg_cx, kegg_norm, KEGG_CMAP, r"KEGG  $-\log_{10}$ FDR")
    cbar(heat_cx, Normalize(-args.vlim, args.vlim), "RdBu_r", "z-score",
         ticks=[-args.vlim, 0, args.vlim])
    cbar(cog_cx, cog_norm, COG_CMAP, r"COG  $-\log_{10}$ FDR")

    for cx, title in [(kegg_cx, "KEGG pathway enrichment"),
                      (heat_cx, "Expression"),
                      (cog_cx, "COG enrichment")]:
        fig.text(cx, title_y, title, ha="center", va="bottom", fontsize=12,
                 fontweight="bold")


def add_bottom_strip(fig, ax_heat_bottom, sample_cols, bounds):
    """Strain (and Chlamy) color strips + rotated sample names, placed BELOW the
    bottom heatmap band (strain strip nearest the heatmap, names at the very
    bottom). Vertical group separators match the heatmap's."""
    hm = ax_heat_bottom.get_position()
    heat_cx = (hm.x0 + hm.x1) / 2
    n_cols = len(sample_cols)
    strains, chlamy, multi_chlamy, strain_pal, chlamy_pal = strain_chlamy_info(sample_cols)

    strip_h = 0.009

    def strip(y, colors):
        ax = fig.add_axes([hm.x0, y, hm.width, strip_h])
        for i, c in enumerate(colors):
            ax.add_patch(plt.Rectangle((i, 0), 1, 1, facecolor=to_rgb(c),
                                       edgecolor="none"))
        ax.set_xlim(0, n_cols); ax.set_ylim(0, 1)
        ax.set_yticks([]); ax.set_xticks([])
        for x in bounds:
            ax.axvline(x, color="black", lw=0.8)
        return ax

    y_strain = hm.y0 - 0.004 - strip_h
    bottom_strip = strip(y_strain, [strain_pal[s] for s in strains])
    if multi_chlamy:
        bottom_strip = strip(y_strain - 0.001 - strip_h,
                             [chlamy_pal[c] for c in chlamy])

    bottom_strip.set_xticks(np.arange(n_cols) + 0.5)
    bottom_strip.set_xticklabels(sample_cols, rotation=90, fontsize=6.5)
    bottom_strip.xaxis.set_ticks_position("bottom")
    bottom_strip.tick_params(axis="x", length=0, pad=1.5)

    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    lows = [inv.transform(lbl.get_window_extent().corners())[:, 1].min()
            for lbl in bottom_strip.get_xticklabels() if lbl.get_text()]
    names_bottom = min(lows) if lows else y_strain
    present = [s for s in STRAIN_ORDER if s in set(strains)]
    leg_h = [Patch(facecolor=strain_pal[s], label=s) for s in present]
    ltitle = "Strain"
    if multi_chlamy:
        leg_h += [Patch(facecolor=chlamy_pal[False], label="Chlamy -"),
                  Patch(facecolor=chlamy_pal[True], label="Chlamy +")]
        ltitle = "Strain / Chlamy"
    fig.legend(handles=leg_h, title=ltitle, loc="upper center",
               bbox_to_anchor=(heat_cx, names_bottom - 0.015), frameon=False,
               fontsize=8, title_fontsize=8, ncol=len(leg_h),
               handletextpad=0.5, columnspacing=1.6)


if __name__ == "__main__":
    main()
