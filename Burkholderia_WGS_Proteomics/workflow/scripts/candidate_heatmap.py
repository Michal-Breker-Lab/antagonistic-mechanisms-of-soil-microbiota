#!/usr/bin/env python3
"""Abundance heatmap of the toxin candidates across strains and conditions."""


import os
import sys
import re
import textwrap
import urllib.parse
import warnings
from types import SimpleNamespace

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
from matplotlib.patches import Patch
from matplotlib.colors import Normalize, to_rgb
from matplotlib.path import Path
import seaborn as sns

NA_COLOR = "black"
MISSING_GREY = "#d0d0d0"
CONTIG_GREY = "#bfbfbf"
CONTIG_HATCH = "///"
OUT_OF_SCALE = 0.72
LFC_CMAP = "RdBu_r"
ABUND_CMAP = "Spectral_r"
STRAIN_ORDER = ["MF6", "27D6", "34F7"]
HDR_PT, HDR_GAP_PT = 22.2, 4.75
HDR_PT_T = 52.0
META_COLS = ["Protein", "Gene", "Protein Length", "Combined Total Peptides",
             "Combined Unique Spectral Count", "Protein Probability",
             "Description"]

SYSID_RE = re.compile(r"^(?P<contig>.+?_\d+)_(?P<system>.+)_(?P<num>\d+)$")

COMPARISONS = [
    ("MF6CvsMF6", "MF6+C\nvs MF6"),
    ("27D6vsMF6", "27D6\nvs MF6"),
    ("34F7vsMF6", "34F7\nvs MF6"),
    ("27D6CvsMF6C", "27D6+C\nvs MF6+C"),
    ("34F7CvsMF6C", "34F7+C\nvs MF6+C"),
]

DE_PATHS = {
    "MF6CvsMF6":   "DE/MF6/DE_coculture_d2.tsv",
    "27D6vsMF6":   "DE/27D6/DE_27D6_alone_vs_MF6_alone.tsv",
    "34F7vsMF6":   "DE/34F7/DE_34F7_alone_vs_MF6_alone.tsv",
    "27D6CvsMF6C": "DE/27D6/DE_27D6_withC_vs_MF6_withC.tsv",
    "34F7CvsMF6C": "DE/34F7/DE_34F7_withC_vs_MF6_withC.tsv",
}

ONOFF_SPEC = {
    "MF6CvsMF6":   ("DE/MF6/on_off_withC_d2_vs_alone_d2.tsv",
                    "MF6_withC_d2", "MF6_alone_d2"),
    "27D6vsMF6":   ("DE/27D6/on_off_27D6_alone_vs_MF6_alone.tsv",
                    "27D6_alone_d2", "MF6_alone_d2"),
    "34F7vsMF6":   ("DE/34F7/on_off_34F7_alone_vs_MF6_alone.tsv",
                    "34F7_alone_d2", "MF6_alone_d2"),
    "27D6CvsMF6C": ("DE/27D6/on_off_27D6_withC_vs_MF6_withC.tsv",
                    "27D6_withC_d2", "MF6_withC_d2"),
    "34F7CvsMF6C": ("DE/34F7/on_off_34F7_withC_vs_MF6_withC.tsv",
                    "34F7_withC_d2", "MF6_withC_d2"),
}

MUTANT_COMPS = {"27D6vsMF6", "34F7vsMF6", "27D6CvsMF6C", "34F7CvsMF6C"}
MUTANT_STRAINS = {"27D6", "34F7"}
ABSENT_IN_MUTANTS = set()


UP_COL = "MF6C_vs_MF6_d2_LFC"
SORT_COMP = "MF6CvsMF6"


def parse_sysid(value):
    """contig_1_T6SSi_19 -> ('T6SSi', '19'); anything else -> (value, '')."""
    m = SYSID_RE.match(value or "")
    return (m.group("system"), m.group("num")) if m else (value or "", "")


def dedupe(ids, what=""):
    seen, uniq, dropped = set(), [], []
    for i in ids:
        (uniq.append(i), seen.add(i)) if i not in seen else dropped.append(i)
    if dropped:
        print(f"note:{what} {len(dropped)} duplicate ID(s) collapsed to first occurrence: "
              f"{', '.join(dict.fromkeys(dropped))}")
    return uniq


def read_splits(args):
    """-> [(split_name|None, ids, {id: label}, [(block_name, [ids])])], caller's order kept."""
    if args.ids:
        return [(None, dedupe(list(args.ids)), {}, [])]

    if args.ids_file:
        ids, labels = [], {}
        with open(args.ids_file) as fh:
            for line in fh:
                line = line.split("#")[0].rstrip()
                if not line.strip():
                    continue
                parts = [c.strip() for c in line.split("\t")]
                ids.append(parts[0])
                if len(parts) > 1 and parts[1]:
                    labels[parts[0]] = parts[1]
        return [(None, dedupe(ids), labels, [])]

    if not args.from_tsv:
        raise SystemExit("give one of --ids / --ids-file / --from-tsv")

    filters = []
    for w in args.where:
        if "=" not in w:
            raise SystemExit(f"--where needs COL=VALUE, got {w!r}")
        k, v = w.split("=", 1)
        filters.append((k, v))
    with open(args.from_tsv) as fh:
        lines = [ln.rstrip("\n") for ln in fh
                 if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        raise SystemExit(f"{args.from_tsv} has no data rows")
    cols = lines[0].split("\t")
    rows = [dict(zip(cols, ln.split("\t"))) for ln in lines[1:]]
    for col in (args.id_col, args.label_col, args.group_col, args.split_sysid):
        if col and col not in rows[0]:
            raise SystemExit(f"column {col!r} not in {args.from_tsv}; "
                             f"columns: {', '.join(rows[0])}")

    block_col = args.split_sysid or args.group_col
    splits = {}
    if not args.no_filter:
        if UP_COL not in rows[0]:
            raise SystemExit(f"{args.from_tsv} has no column {UP_COL!r}; "
                             f"pass --no-filter to plot every row")
        before = len(rows)
        rows = [r for r in rows
                if is_up_in_coculture(r.get(UP_COL), args.onoff_tier)]
        print(f"rows: {len(rows)} of {before} are ON or significantly up in {UP_COL}"
              + (f"  (ON calls gated at tier >= {args.onoff_tier})"
                 if args.onoff_tier != "weak" else ""))
    for r in rows:
        if any(r.get(k) != v for k, v in filters):
            continue
        pid = (r.get(args.id_col) or "").strip()
        if not pid:
            continue
        sysname, sysnum = ("", "")
        if args.split_sysid:
            sysname, sysnum = parse_sysid((r.get(args.split_sysid) or "").strip())
        key = sysname if args.split_sysid else None
        ids, labels, blocks = splits.setdefault(key, ([], {}, []))
        ids.append(pid)
        if args.label_col and r.get(args.label_col):
            lab = r[args.label_col].strip()
            labels[pid] = f"{lab} [{sysnum}]" if sysnum else lab
        if block_col:
            g = (r.get(block_col) or "").strip()
            if not blocks or blocks[-1][0] != g:
                blocks.append((g, []))
            blocks[-1][1].append(pid)

    if not splits:
        raise SystemExit("the requested list is empty (check --where / --id-col)")
    out = []
    for name, (ids, labels, blocks) in splits.items():
        out.append((name, dedupe(ids, f" [{name}]" if name else ""), labels, blocks))
    return out


def load_de(root, comp, sig_col):
    """This project's limma tables, renamed to the columns the panels expect:
    logFC -> diff, adj.P.Val -> adj_pval / P.Value -> pval, Gene -> Preferred_name."""
    path = os.path.join(root, DE_PATHS[comp])
    if not os.path.exists(path):
        raise SystemExit(f"missing DE table: {path}")
    d = pd.read_csv(path, sep="\t", comment="#").set_index("Protein")
    d = d.rename(columns={"logFC": "diff", "adj.P.Val": "adj_pval",
                          "P.Value": "pval", "Gene": "Preferred_name"})
    if sig_col not in d.columns:
        raise SystemExit(f"{path} has no column for {sig_col!r}")
    return d


def load_onoff(root, comp):
    """-> {protein: 'ON'|'OFF'} for one contrast, oriented like its logFC: ON = absent
    from the denominator group and present in the numerator (so it belongs at the red
    end), OFF = the reverse."""
    path_rel, ga, gb = ONOFF_SPEC[comp]
    path = os.path.join(root, path_rel)
    if not os.path.exists(path):
        print(f"note: no ON/OFF table for {comp} ({path_rel}), those cells stay black")
        return {}
    out = {}
    for r in pd.read_csv(path, sep="\t", comment="#", dtype=str).to_dict("records"):
        pres = (r.get("present_in") or "").strip()
        full = pres if pres.startswith(("MF6_", "27D6_", "34F7_")) else f"MF6_{pres}"
        state = "ON" if full == ga else "OFF" if full == gb else None
        if state:
            out[r["Protein"]] = state
    return out


def contig_of(root, gff):
    """{locus_tag: contig} from the CDS lines of the annotation GFF."""
    path = os.path.join(root, gff)
    if not os.path.exists(path):
        print(f"note: {gff} not found, contig_2 cells cannot be marked")
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"locus_tag=([^;]+)", f[8])
            if m:
                out[m.group(1)] = f[0]
    return out


def sample_groups(cols):
    """Ordered (strain, chlamy) -> [sample columns], strain-major, -Chl before +Chl.
    This project names samples <strain>_<alone|withC>_d2_<rep>."""
    groups = {}
    for c in cols:
        groups.setdefault((c.split("_")[0], "_withC_" in c), []).append(c)
    order = [(s, ch) for s in STRAIN_ORDER for ch in (False, True) if (s, ch) in groups]
    return order, groups


ONOFF_TIERS = ("weak", "moderate", "strong")


def is_up_in_coculture(cell, min_tier="weak"):
    """True for the rows we plot: an ON call, or a significant positive logFC of
    any magnitude.  The cell is candidate_list.tsv's merged contrast cell, so it
    is one of '2.579*' / '-0.099' / 'ON, strong (log2 LFQ 24.74)' / 'OFF, ...' /
    'ND (n/n vs n/n)' / 'NA (contig 2)'.

    min_tier gates the ON calls only - a significant logFC is evidence in its own
    right and is never filtered by it.  An ON call at the bottom of the other group's
    dynamic range may mean "below detection" rather than "switched off", so requiring
    'moderate' (pct > 20) or 'strong' (pct > 70) drops the calls that cannot be
    distinguished from a detection failure.  See tests/detection_probability.md."""
    v = (cell or "").strip()
    if v.startswith("ON,") or v == "ON":
        floor = ONOFF_TIERS.index(min_tier)
        tier = next((x for x in ONOFF_TIERS if f", {x}" in v), None)
        if tier is None:
            return floor == 0
        return ONOFF_TIERS.index(tier) >= floor
    if not v.endswith("*"):
        return False
    try:
        return float(v[:-1]) > 0
    except ValueError:
        return False


def clean(v):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    s = str(v).strip()
    return "" if s in ("NA", "nan", "-", "None") else s


def finite_or(value, fallback):
    return fallback if value is None or not np.isfinite(value) else float(value)


def nanpct(a, q):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanpercentile(np.abs(a), q) if np.isfinite(a).any() else np.nan


ANNOT_WIDTH = 1.0
ANNOT_CELL, ANNOT_HDR = "#f5f5f5", "#dcdcdc"
PRESENT_COLOR = "#1b7837"

CHECK = Path([(-1.0, 0.05), (-0.3, -0.75), (1.0, 0.8)],
             [Path.MOVETO, Path.LINETO, Path.LINETO])


def col_width_pt(ax):
    """Width of one heatmap column, in points."""
    (x0, _), (x1, _) = ax.transData.transform([(0, 0), (1, 0)])
    return abs(x1 - x0) / ax.figure.dpi * 72.0


def cat_cell_pt(ax, args):
    """Extent of ONE category cell along the axis the categories run on, in points -
    the y axis once transposed, the x axis otherwise.  What an annotation label has to
    fit into across its short side."""
    a, b = ((0, 0), (0, 1)) if args.transpose else ((0, 0), (1, 0))
    (x0, y0), (x1, y1) = ax.transData.transform([a, b])
    return abs((y1 - y0) if args.transpose else (x1 - x0)) / ax.figure.dpi * 72.0


def fit_label(text, width_pt, height_pt=None, max_fs=8.5, min_fs=4.5):
    """Wrap `text` and shrink it until it fits a box this wide (and, if given, this
    tall).  The annotation column is only one column wide, so a long label has to be
    laid out to that width rather than assumed to fit."""
    best, fs = (text, min_fs), max_fs
    while fs >= min_fs:
        cap = max(1, int(width_pt / (fs * 0.52)))
        lines = textwrap.wrap(text, width=cap, break_long_words=False) or [text]
        wide = max(len(l) for l in lines) * fs * 0.52
        tall = len(lines) * fs * 1.18
        if wide <= width_pt and (height_pt is None or tall <= height_pt):
            return "\n".join(lines), fs
        best = ("\n".join(lines), fs)
        fs -= 0.25
    return best


def hdr_pt(args):
    """Thickness of the header band in points (a height, or a width once transposed)."""
    return HDR_PT_T if args.transpose else HDR_PT


def hdr_geom(args):
    """(offset, thickness, text centre) of the header band, in cell units.  Cells are
    args.row_height inches along whichever axis the band crosses, in both layouts."""
    rows_per_pt = 1.0 / (args.row_height * 72.0)
    h = hdr_pt(args) * rows_per_pt
    y = -0.5 - HDR_GAP_PT * rows_per_pt - h
    return y, h, y + h / 2


def annot_x(n_cols):
    """Centre of the annotation column."""
    return n_cols - 0.5 + ANNOT_WIDTH / 2


def annot_right(n_cols):
    return n_cols - 0.5 + ANNOT_WIDTH


def cell_xy(j, i, T):
    """(x, y) of the cell in category j, protein i.  Transposing the figure is exactly
    this swap - proteins run along x instead of y - so every draw call goes through it
    and no panel can end up half-transposed."""
    return (i, j) if T else (j, i)


def new_axes(n_rows, n_cols, labels, catlabels, args, width_per_col, pad, annot=None):
    """n_rows = proteins, n_cols = contrasts / sample groups.  `labels` names the
    proteins, `catlabels` the categories; which axis each lands on is --transpose."""
    T = args.transpose
    far = annot_right(n_cols) if annot else n_cols - 0.5
    if T:
        fig, ax = plt.subplots(figsize=(args.row_height * n_rows + 3.4,
                                        args.cat_height * (far + 0.5) + 3.0))
        ax.set_xlim(-0.5, n_rows - 0.5)
        ax.set_ylim(far, -0.5)
    else:
        fig, ax = plt.subplots(figsize=(width_per_col * (far + 0.5) + 5.4,
                                        args.row_height * n_rows + pad))
        ax.set_xlim(-0.5, far)
        ax.set_ylim(n_rows - 0.5, -0.5)

    foot = annot and not annot.get("header")
    pos, lab, fs = list(range(n_cols)), list(catlabels), None
    if foot:
        if T:
            txt, fs = "\n".join(textwrap.wrap(annot["label"], 16)), 8
        else:
            txt, fs = fit_label(annot["label"], col_width_pt(ax) * ANNOT_WIDTH)
        pos.append(annot_x(n_cols)); lab.append(txt)

    cat_ax, prot_ax = (ax.yaxis, ax.xaxis) if T else (ax.xaxis, ax.yaxis)
    cat_ax.set_ticks(pos)
    cat_ax.set_ticklabels(lab, fontsize=8)
    if fs:
        cat_ax.get_ticklabels()[-1].set_fontsize(fs)
    prot_ax.set_ticks(range(n_rows))
    prot_ax.set_ticklabels(labels, fontsize=args.label_fontsize,
                           **(dict(rotation=45, ha="left", va="bottom",
                                   rotation_mode="anchor") if T else {}))
    if T:
        prot_ax.set_ticks_position("top")
        prot_ax.set_label_position("top")
        cat_ax.set_ticks_position("right")
        cat_ax.set_label_position("right")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_color("black"); sp.set_linewidth(0.8)
    return fig, ax


def draw_annot(ax, n_cols, n_rows, args, annot):
    """The presence/absence strip: one cell per row, ruled exactly like the heatmap's
    own cells, a green check where the protein is present and blank where it is not.
    Categorical and outside the heatmap's colour scale by construction."""
    if not annot:
        return
    T = args.transpose
    c = annot_x(n_cols)
    ms = float(np.clip(args.row_height * 72.0 * 0.60, 5.0, 16.0))
    (ax.axhline if T else ax.axvline)(n_cols - 0.5, color="black", lw=1.0, zorder=3)
    for i in range(n_rows):
        xy = (i - 0.5, n_cols - 0.5) if T else (n_cols - 0.5, i - 0.5)
        w, h = (1, ANNOT_WIDTH) if T else (ANNOT_WIDTH, 1)
        ax.add_patch(plt.Rectangle(xy, w, h, facecolor=ANNOT_CELL, edgecolor="black",
                                   linewidth=0.5, zorder=1))
        if annot["values"][i]:
            x, y = cell_xy(c, i, T)
            ax.plot(x, y, marker=CHECK, markersize=ms, linestyle="none",
                    markerfacecolor="none", markeredgecolor=annot["color"],
                    markeredgewidth=max(1.2, ms * 0.16), zorder=4)
    if annot.get("header"):
        hy, hh, hty = hdr_geom(args)
        xy = (hy, n_cols - 0.5) if T else (n_cols - 0.5, hy)
        w, h = (hh, ANNOT_WIDTH) if T else (ANNOT_WIDTH, hh)
        ax.add_patch(plt.Rectangle(xy, w, h, facecolor=ANNOT_HDR, edgecolor="black",
                                   lw=0.8, clip_on=False, zorder=0))
        box = ((hdr_pt(args), cat_cell_pt(ax, args) * ANNOT_WIDTH) if T else
               (col_width_pt(ax) * ANNOT_WIDTH, HDR_PT))
        txt, fs = fit_label(annot["label"], *box)
        tx, ty = cell_xy(c, hty, T)
        ax.text(tx, ty, txt, ha="center", va="center", fontsize=fs, linespacing=1.1,
                clip_on=False, zorder=1)


def draw_blocks(ax, blocks, n_cols, args, annot=None):
    T = args.transpose
    far = (ax.get_ylim()[0] if T else ax.get_xlim()[1]) - (0.08 if annot else 0.42)
    for gname, r0, r1 in blocks:
        if r0 > 0:
            (ax.axvline if T else ax.axhline)(r0 - 0.5, color="black", lw=1.2, zorder=3)
        x, y = cell_xy(far, (r0 + r1) / 2, T)
        ax.text(x, y, gname, rotation=0 if T else 270, va="top" if T else "center",
                ha="center" if T else "left", fontsize=7, clip_on=False)


def cbar_pad(ax, args, catlabels, n_x, blocks):
    """Gap between the matrix and the colourbar, as a fraction of the axes width.
    Transposed, the category labels move to the right edge and the colourbar has to
    clear them, so the gap is measured from the widest label rather than fixed."""
    if not args.transpose:
        return 0.10 if blocks else 0.02
    wid = max((max(len(l) for l in str(lab).split("\n")) for lab in catlabels), default=4)
    need_pt = wid * 8 * 0.52 + 12
    return max(0.02, need_pt / max(col_width_pt(ax) * n_x, 1.0))


def place_legend(ax, handles, labels, args, ncol=2, above_frac=0.0):
    """Above the matrix normally; BELOW it when transposed, where the top edge is taken
    by the rotated protein names.  `above_frac` clears the strain strip in the untransposed
    abundance panel."""
    if args.transpose:
        return ax.legend(handles=handles, labels=labels, loc="upper left",
                         bbox_to_anchor=(0.0, -0.04), ncol=ncol, frameon=False,
                         fontsize=8, handlelength=1.6, borderpad=0.4)
    return ax.legend(handles=handles, labels=labels, loc="lower left",
                     bbox_to_anchor=(0.0, 1.001 + above_frac), ncol=ncol, frameon=False,
                     fontsize=8, handlelength=1.6, borderpad=0.4)


def cap_line(args):
    """The header note for a clamped abundance colour scale.  The MEANS IN THIS TABLE
    ARE NEVER CLAMPED - the limit is a drawing choice, so it has to be stated here or a
    reader comparing table to figure would think the numbers had been altered."""
    lo, hi = args.abund_vmin, args.abund_vmax
    if lo is None and hi is None:
        return ""
    what = (f"capped at {hi:g}" if lo is None else
            f"floored at {lo:g}" if hi is None else
            f"clamped to {lo:g}-{hi:g}")
    return (f"# figure colour scale {what} (--abund-vmin/--abund-vmax); the means below "
            "are the real values, unclamped - cells past the limit are drawn in the end "
            "colour and the colourbar carries an extend cap\n")


def vector_colorbar(cbar):
    """matplotlib rasterizes the colourbar's gradient by default, which lands in the PDF
    as a ~100 ppi image among otherwise pure vector art.  Turn it off so the whole
    figure scales cleanly and Illustrator gets shapes rather than a bitmap."""
    if getattr(cbar, "solids", None) is not None:
        cbar.solids.set_rasterized(False)
    for patch in getattr(cbar, "solids_patches", []) or []:
        patch.set_rasterized(False)


def save(fig, args, prefix, suffix):
    for ext in args.formats:
        out = os.path.join(args.outdir, f"{prefix}_{suffix}.{ext}")
        fig.savefig(out, bbox_inches="tight")
        print(f"  wrote {out}")
    if args.png:
        png = os.path.join(args.outdir, f"{prefix}_{suffix}.png")
        fig.savefig(png, bbox_inches="tight", dpi=150)
        print(f"  wrote {png}")
    plt.close(fig)


def plot_lfc(diffmat, sig_cell, onoff, absent, labels, blocks, args, prefix, title, vlim,
             annot=None):
    n_rows, n_cols = diffmat.shape
    norm = Normalize(vmin=-vlim, vmax=vlim)
    plt.rcParams["hatch.linewidth"] = 0.4

    cmap = plt.get_cmap(args.cmap).copy()
    cmap.set_over([c * OUT_OF_SCALE for c in to_rgb(cmap(1.0))])
    cmap.set_under([c * OUT_OF_SCALE for c in to_rgb(cmap(0.0))])
    disp = diffmat.copy()
    disp[onoff == "ON"] = vlim * 1.2
    disp[onoff == "OFF"] = -vlim * 1.2
    is_lfc = np.isfinite(diffmat)
    ms = float(np.clip(args.row_height * 72.0 * 0.35, 3.0, 7.0))

    T = args.transpose
    fig, ax = new_axes(n_rows, n_cols, labels, [lab for _, lab in COMPARISONS],
                       args, 0.85, 2.2, annot)
    for i in range(n_rows):
        for j in range(n_cols):
            x, y = cell_xy(j, i, T)
            if absent[i, j]:
                ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                           facecolor=CONTIG_GREY, edgecolor="black",
                                           linewidth=0.5, hatch=CONTIG_HATCH, zorder=1))
                continue
            na = np.isnan(disp[i, j])
            color = NA_COLOR if na else cmap(norm(disp[i, j]))
            hatch = None if (not is_lfc[i, j] or sig_cell[i, j]) else "xxx"
            ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=color,
                                       edgecolor="black", linewidth=0.5,
                                       hatch=hatch, zorder=1))
            if onoff[i, j]:
                ax.plot(x, y, marker="^" if onoff[i, j] == "ON" else "v", markersize=ms,
                        markerfacecolor="white", markeredgecolor="white",
                        markeredgewidth=0.4, linestyle="none", zorder=4)
            if args.annotate and is_lfc[i, j]:
                ax.text(x, y, f"{diffmat[i, j]:.1f}", ha="center", va="center",
                        fontsize=5.5, zorder=4,
                        color="white" if abs(diffmat[i, j]) > 0.65 * vlim else "black")

    for v in (0.5, 2.5):
        (ax.axhline if T else ax.axvline)(v, color="black", lw=1.2, zorder=3)
    draw_annot(ax, n_cols, n_rows, args, annot)
    draw_blocks(ax, blocks, n_cols, args, annot)

    has_on, has_off = (onoff == "ON").any(), (onoff == "OFF").any()
    extend = ("both" if has_on and has_off else
              "max" if has_on else "min" if has_off else "neither")
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                        fraction=0.045, extend=extend,
                        pad=cbar_pad(ax, args, [lab for _, lab in COMPARISONS]
                                     + ([annot["label"]] if annot else []),
                                     n_rows, blocks))
    vector_colorbar(cbar)
    cbar.set_label(r"$\log_2$ FC (strain/condition $-$ MF6)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    for present, y, va, mark, name in ((has_on, 1.0, "bottom", "^", "ON"),
                                       (has_off, 0.0, "top", "v", "OFF")):
        if present:
            cbar.ax.plot(0.5, y + (0.030 if va == "bottom" else -0.030), marker=mark,
                         markersize=ms, markerfacecolor="white", markeredgecolor="white",
                         markeredgewidth=0.4, linestyle="none", clip_on=False,
                         transform=cbar.ax.transAxes, zorder=5)
            cbar.ax.text(0.5, y + (0.075 if va == "bottom" else -0.075), name,
                         transform=cbar.ax.transAxes, ha="center", va=va, fontsize=8)

    handles = [Patch(facecolor="white", edgecolor="black", hatch="xxx")]
    labels = [f"not significant ({args.sig_col} ≥ {args.sig:g})"]
    if absent.any():
        handles.append(Patch(facecolor=CONTIG_GREY, edgecolor="black", hatch=CONTIG_HATCH))
        labels.append(f"{'/'.join(sorted(ABSENT_IN_MUTANTS))} absent in mutant")
    handles.append(Patch(facecolor=NA_COLOR, edgecolor="black"))
    labels.append("NA (not quantified)")
    place_legend(ax, handles, labels, args, ncol=min(len(handles), 3))
    if title:
        ax.set_title(title, fontsize=10, pad=26)
    save(fig, args, prefix, "lfc_heatmap")


def plot_abundance(mean, in_matrix, absent_group, order, labels, blocks, args, prefix,
                   title, annot=None):
    annot = dict(annot, header=True) if annot else None
    n_rows, n_cols = mean.shape
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        dmin = np.nanmin(mean) if np.isfinite(mean).any() else 0.0
        dmax = np.nanmax(mean) if np.isfinite(mean).any() else 1.0
    vmin = dmin if args.abund_vmin is None else args.abund_vmin
    vmax = dmax if args.abund_vmax is None else args.abund_vmax
    if vmin == vmax:
        vmin, vmax = vmin - 0.5, vmax + 0.5
    lo_cap = args.abund_vmin is not None and dmin < vmin
    hi_cap = args.abund_vmax is not None and dmax > vmax
    abund_extend = ("both" if lo_cap and hi_cap else
                    "max" if hi_cap else "min" if lo_cap else "neither")
    norm, cmap = Normalize(vmin=vmin, vmax=vmax), plt.get_cmap(args.abund_cmap)

    T = args.transpose
    fig, ax = new_axes(n_rows, n_cols, labels,
                       ["Chl−" if not ch else "Chl+" for _, ch in order],
                       args, 0.55, 2.0, annot)
    plt.rcParams["hatch.linewidth"] = 0.4
    for i in range(n_rows):
        for j in range(n_cols):
            x, y = cell_xy(j, i, T)
            if absent_group[i, j]:
                ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                           facecolor=CONTIG_GREY, edgecolor="black",
                                           linewidth=0.5, hatch=CONTIG_HATCH, zorder=1))
                continue
            if not in_matrix[i]:
                color = NA_COLOR
            elif np.isnan(mean[i, j]):
                color = MISSING_GREY
            else:
                color = cmap(norm(mean[i, j]))
            ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor=color,
                                       edgecolor="black", linewidth=0.5, zorder=1))
            if args.annotate and in_matrix[i] and not np.isnan(mean[i, j]):
                ax.text(x, y, f"{mean[i, j]:.0f}", ha="center", va="center",
                        fontsize=5.5, zorder=4, color="black")

    hy, hh, hty = hdr_geom(args)
    strain_pal = dict(zip(STRAIN_ORDER, sns.color_palette("Set2", 3)))
    for j in range(1, n_cols):
        if order[j][0] != order[j - 1][0]:
            (ax.axhline if T else ax.axvline)(j - 0.5, color="black", lw=1.2, zorder=3)
    j = 0
    while j < n_cols:
        s, k = order[j][0], j
        while k < n_cols and order[k][0] == s:
            k += 1
        xy = (hy, j - 0.5) if T else (j - 0.5, hy)
        w, h = (hh, k - j) if T else (k - j, hh)
        ax.add_patch(plt.Rectangle(xy, w, h, clip_on=False,
                                   facecolor=strain_pal[s], edgecolor="black", lw=0.8))
        tx, ty = cell_xy((j + k - 1) / 2, hty, T)
        ax.text(tx, ty, s, ha="center", va="center",
                fontsize=10, fontweight="bold", clip_on=False)
        j = k
    draw_annot(ax, n_cols, n_rows, args, annot)
    draw_blocks(ax, blocks, n_cols, args, annot)

    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                        fraction=0.045, extend=abund_extend,
                        pad=cbar_pad(ax, args,
                                     ["Chl−" if not ch else "Chl+" for _, ch in order],
                                     n_rows, blocks))
    vector_colorbar(cbar)
    cbar.set_label(r"mean $\log_2$ normalized LFQ", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    handles = [Patch(facecolor=MISSING_GREY, edgecolor="black",
                     label="not detected in that group")]
    if absent_group.any():
        handles.append(Patch(facecolor=CONTIG_GREY, edgecolor="black", hatch=CONTIG_HATCH,
                             label=f"{'/'.join(sorted(ABSENT_IN_MUTANTS))} absent in mutant"))
    if not in_matrix.all():
        handles.append(Patch(facecolor=NA_COLOR, edgecolor="black",
                             label="not in the pellet matrix"))
    ax_h_pt = max(args.row_height * n_rows * 72.0, 1.0)
    strip_pt = 0.0 if T else abs(hy) * args.row_height * 72.0
    place_legend(ax, handles, [h.get_label() for h in handles], args,
                 above_frac=(strip_pt + 4) / ax_h_pt)
    if title:
        ax.set_title(title, fontsize=10, pad=strip_pt + 30)
    save(fig, args, prefix, "abundance_heatmap")


def sort_key(B, kind):
    """The quantity a panel is ordered on, larger = higher up.

    lfc        the SORT_COMP logFC, with a presence/absence call taken to its limit:
               ON = +inf (above every finite fold change), OFF = -inf (below every
               one).  NaN for a cell with no call at all, which sorts last.
    abundance  mean of the group means, i.e. how much protein there is overall.
    """
    if kind == "abundance":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(B["mean"], axis=1)
    c = [i for i, (comp, _) in enumerate(COMPARISONS) if comp == SORT_COMP]
    if not c:
        raise SystemExit(f"--sort lfc needs {SORT_COMP!r} among the contrasts")
    k = B["diffmat"][:, c[0]].astype(float).copy()
    k[B["onoff"][:, c[0]] == "ON"] = np.inf
    k[B["onoff"][:, c[0]] == "OFF"] = -np.inf
    return k


def row_order(B, key, contig_last):
    """Permutation for one panel: blocks stay contiguous and in place (--group-col /
    --split-sysid still work), the mutants' missing replicon is held in a block at the
    bottom of each, then `key` descending with no-value rows last, stable on ties."""
    blk = B["block_idx"]
    tail = B["on_contig"].astype(int) if contig_last else np.zeros(len(key), int)
    return sorted(range(len(key)),
                  key=lambda i: (blk[i], tail[i], int(np.isnan(key[i])),
                                 0.0 if np.isnan(key[i]) else -key[i], i))


def sort_kind(args, panel):
    """Which key this panel sorts on, given --sort."""
    return panel if args.sort == "panel" else args.sort


def sort_note(args, panel):
    """One line for the TSV header saying how the rows ended up in this order."""
    kind = sort_kind(args, panel)
    what = {"lfc": f"descending {SORT_COMP} logFC (ON first, then the fold changes, "
                   f"then OFF, then rows with no call)",
            "abundance": "descending overall intensity (mean of the group means)",
            "input": "the order given by the input list (not sorted)"}[kind]
    tail = ("" if args.no_contig_last or kind == "input" else
            f"; {'/'.join(sorted(ABSENT_IN_MUTANTS))} genes held in a block at the "
            "bottom, sorted among themselves by the same key")
    return f"# rows: {what}{tail}\n"


def panel_view(B, args, panel):
    """B with its rows in the order this panel should draw them."""
    kind = sort_kind(args, panel)
    if kind == "input" and args.no_contig_last:
        return B
    key = (np.zeros(len(B["ids"])) if kind == "input" else sort_key(B, kind))
    return reorder(B, row_order(B, key, not args.no_contig_last))


def reorder(B, perm):
    """A copy of the built panel data with its rows permuted.  Blocks are recomputed
    from the permuted block indices - they stay contiguous because the block index is
    row_order's most significant key."""
    out = dict(B)
    for k in ("diffmat", "pmat", "sig_cell", "onoff", "absent", "mean", "absent_group",
              "in_matrix", "on_contig", "block_idx"):
        out[k] = B[k][perm]
    for k in ("ids", "labels"):
        out[k] = [B[k][i] for i in perm]
    blocks, bi = [], out["block_idx"]
    for r, b in enumerate(bi):
        if blocks and bi[blocks[-1][1]] == b:
            blocks[-1] = (blocks[-1][0], blocks[-1][1], r)
        else:
            blocks.append((B["block_names"].get(b, ""), r, r))
    out["blocks"] = [x for x in blocks if x[0]]
    return out


def build(ids, user_labels, raw_blocks, args, de, onoff_tbl, contig, mat, sgroup_order,
          sgroups, want_lfc, want_abund):
    """Assemble every matrix/label the two panels need for one split."""
    n_rows, n_cols = len(ids), len(COMPARISONS)
    diffmat = np.full((n_rows, n_cols), np.nan)
    pmat = np.full((n_rows, n_cols), np.nan)
    sig_cell = np.zeros((n_rows, n_cols), dtype=bool)
    onoff = np.full((n_rows, n_cols), "", dtype=object)
    absent = np.zeros((n_rows, n_cols), dtype=bool)
    for c, (comp, _) in enumerate(COMPARISONS):
        d = de[comp]
        mutant = comp in MUTANT_COMPS
        for i, pid in enumerate(ids):
            if mutant and contig.get(pid) in ABSENT_IN_MUTANTS:
                absent[i, c] = True
                continue
            onoff[i, c] = onoff_tbl[comp].get(pid, "")
            if pid not in d.index:
                continue
            diffmat[i, c] = d.at[pid, "diff"]
            pmat[i, c] = d.at[pid, args.sig_col]
            sig_cell[i, c] = pd.notna(pmat[i, c]) and pmat[i, c] < args.sig

    in_matrix = np.array([pid in mat.index for pid in ids], dtype=bool)
    mean = np.full((n_rows, len(sgroup_order)), np.nan)
    for j, g in enumerate(sgroup_order):
        vals = np.full((n_rows, len(sgroups[g])), np.nan)
        for i, pid in enumerate(ids):
            if in_matrix[i]:
                vals[i] = mat.loc[pid, sgroups[g]].to_numpy(float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            mean[:, j] = np.nanmean(vals, axis=1)

    absent_group = np.zeros((n_rows, len(sgroup_order)), dtype=bool)
    for j, (s, _ch) in enumerate(sgroup_order):
        if s in MUTANT_STRAINS:
            for i, pid in enumerate(ids):
                absent_group[i, j] = contig.get(pid) in ABSENT_IN_MUTANTS
    mean[absent_group] = np.nan

    gene, product = {}, {}
    for comp, _ in COMPARISONS:
        d = de[comp]
        for pid in ids:
            if pid in d.index and pid not in gene:
                gene[pid] = clean(d.at[pid, "Preferred_name"])
                product[pid] = clean(d.at[pid, "Description"])
    for i, pid in enumerate(ids):
        if pid not in gene and in_matrix[i]:
            gene[pid] = clean(mat.at[pid, "Gene"])
            product[pid] = clean(mat.at[pid, "Description"])

    has_lfc = ~np.isnan(diffmat).all(axis=1) | (onoff != "").any(axis=1)
    has_abund = ~np.isnan(mean).all(axis=1) | absent_group.any(axis=1)
    stats = (int(has_lfc.sum()), int(has_abund.sum()))
    if args.drop_all_na:
        keep = (has_lfc if want_lfc else np.zeros(n_rows, bool)) | \
               (has_abund if want_abund else np.zeros(n_rows, bool))
        ids = [p for p, k in zip(ids, keep) if k]
        diffmat, pmat, sig_cell = diffmat[keep], pmat[keep], sig_cell[keep]
        onoff, absent = onoff[keep], absent[keep]
        mean, in_matrix, absent_group = mean[keep], in_matrix[keep], absent_group[keep]
        n_rows = len(ids)

    labels = []
    for pid in ids:
        name = user_labels.get(pid) or gene.get(pid) or \
            urllib.parse.unquote(product.get(pid, "")) or pid
        full = f"{name} ({str(pid).split('_')[-1]})"
        labels.append(full if args.transpose else
                      "\n".join(textwrap.wrap(full, width=40, break_long_words=False,
                                              break_on_hyphens=False) or [full]))

    pos = {p: i for i, p in enumerate(ids)}
    blocks = []
    for gname, gids in raw_blocks:
        rows = sorted(pos[p] for p in dict.fromkeys(gids) if p in pos)
        if rows:
            blocks.append((gname, rows[0], rows[-1]))

    gname_of = {p: gname for gname, gids in raw_blocks for p in gids}
    block_idx = np.zeros(n_rows, dtype=int)
    block_names, run, prev = {}, 0, None
    for i, pid in enumerate(ids):
        g = gname_of.get(pid, "")
        if i and g != prev:
            run += 1
        block_idx[i], block_names[run], prev = run, g, g

    on_contig = np.array([contig.get(p) in ABSENT_IN_MUTANTS for p in ids], dtype=bool)

    return dict(ids=ids, diffmat=diffmat, pmat=pmat, sig_cell=sig_cell, onoff=onoff,
                absent=absent, mean=mean, absent_group=absent_group,
                in_matrix=in_matrix, on_contig=on_contig, block_idx=block_idx,
                block_names=block_names, gene=gene, product=product, labels=labels,
                blocks=blocks, stats=stats)


def write_tsv(B, args, prefix, user_labels, order, kind, src_line, annot=None):
    ids, blocks, n_rows = B["ids"], B["blocks"], len(B["ids"])
    t = pd.DataFrame({
        "protein_id": ids,
        "label": [user_labels.get(p, "") for p in ids],
        "gene": [B["gene"].get(p, "") for p in ids],
        "product": [urllib.parse.unquote(B["product"].get(p, "")) for p in ids],
    })
    if annot:
        t[annot["column"]] = ["yes" if v else "no" for v in annot["values"]]
        src_line = src_line + annot["note"]
    if blocks:
        gof = {i: g for g, r0, r1 in blocks for i in range(r0, r1 + 1)}
        t.insert(1, "group", [gof.get(i, "") for i in range(n_rows)])
    if not any(t["label"]):
        t = t.drop(columns="label")
    order_line = sort_note(args, kind)

    if kind == "lfc":
        for c, (comp, _) in enumerate(COMPARISONS):
            t[f"{comp}_logFC"] = B["diffmat"][:, c]
            t[f"{comp}_{args.sig_col}"] = B["pmat"][:, c]
            t[f"{comp}_call"] = [
                "absent (" + "/".join(sorted(ABSENT_IN_MUTANTS)) + ")" if a
                else s or ("LFC" if np.isfinite(v) else "ND")
                for a, s, v in zip(B["absent"][:, c], B["onoff"][:, c], B["diffmat"][:, c])]
        path = os.path.join(args.outdir, f"{prefix}_lfc_heatmap.tsv")
        head = (f"# log2 FC vs MF6 for the {n_rows} proteins in {prefix}_lfc_heatmap.svg\n"
                + src_line +
                "# <contrast>_logFC = limma diff (strain/condition - MF6): >0 higher in the "
                "other, <0 higher in MF6; blank whenever <contrast>_call is not LFC\n"
                f"# <contrast>_{args.sig_col} = limma {args.sig_col}; significant when "
                f"< {args.sig:g} (cells with {args.sig_col} >= {args.sig:g} are cross-hatched)\n"
                "# <contrast>_call = how the cell is drawn: LFC = the number above; "
                "ON = never detected in the reference group but present here (scale-extreme "
                "red + up triangle); OFF = the reverse (scale-extreme blue + down triangle); "
                "absent = replicon deleted in the mutant, no comparison exists (grey hatch); "
                "ND = not testable in that contrast (black).  ON/OFF have no finite fold "
                "change - the fill marks direction only, never magnitude\n"
                + order_line)
    else:
        t["in_pellet_matrix"] = ["yes" if v else "no" for v in B["in_matrix"]]
        for j, (s, ch) in enumerate(order):
            t[f"{s}_{'Chl+' if ch else 'Chl-'}_meanLog2LFQ"] = B["mean"][:, j]
        path = os.path.join(args.outdir, f"{prefix}_abundance_heatmap.tsv")
        head = (f"# mean log2 normalized LFQ for the {n_rows} proteins in "
                f"{prefix}_abundance_heatmap.svg\n" + src_line +
                f"# source: {args.data} (day {args.day}"
                + (f", excluded: {', '.join(args.exclude)}" if args.exclude else "") + ")\n"
                "# <strain>_<Chl+/->_meanLog2LFQ = mean over the group's detected replicates; "
                "blank = detected in no replicate of that group (grey cell)\n"
                f"# {'/'.join(sorted(ABSENT_IN_MUTANTS))} genes are censored (blank, grey "
                "hatch) in the 27D6/34F7 columns: both mutants have lost that replicon, so "
                "any MaxLFQ value there is a match-between-runs transfer from the wild type\n"
                "# in_pellet_matrix = no -> the protein is absent from the matrix entirely "
                "(black row in the figure)\n"
                + cap_line(args) + order_line)
    with open(path, "w") as fh:
        fh.write(head)
        t.to_csv(fh, sep="\t", index=False)
    print(f"  wrote {path}")


def main(args):
    """Draw the candidate heatmap panels and write their TSVs.

    The panel drawing still accepts a categorical row annotation (`draw_annot`);
    nothing supplies one since the MF6-supernatant column was dropped, so `a` is
    always None here."""
    global ABSENT_IN_MUTANTS
    ABSENT_IN_MUTANTS = {args.absent_replicon}
    splits = read_splits(args)
    os.makedirs(args.outdir, exist_ok=True)
    want_lfc = args.mode in ("both", "lfc")
    want_abund = args.mode in ("both", "abundance")

    de = {comp: load_de(args.root, comp, args.sig_col) for comp, _ in COMPARISONS}
    onoff_tbl = {comp: ({} if args.no_onoff else load_onoff(args.root, comp))
                 for comp, _ in COMPARISONS}
    contig = contig_of(args.root, args.gff)
    mat = pd.read_csv(os.path.join(args.root, args.data), sep="\t").set_index("Protein")
    keep_cols = [c for c in mat.columns if c not in META_COLS
                 and re.search(r"_d2_[1-4]$", c) and c not in set(args.exclude)]
    if want_abund and not keep_cols:
        raise SystemExit(f"no day-{args.day} sample columns left after exclusions")
    order, sgroups = sample_groups(keep_cols)

    built = []
    for name, ids, labels, raw_blocks in splits:
        B = build(ids, labels, raw_blocks, args, de, onoff_tbl, contig, mat, order,
                  sgroups, want_lfc, want_abund)
        built.append((name, labels, B))

    shared = None
    if args.shared_vlim:
        allv = np.concatenate([B["diffmat"].ravel() for _, _, B in built])
        shared = finite_or(nanpct(allv, 98), 1.0)
        print(f"shared LFC colour limit: ±{shared:.2f}")

    src = args.from_tsv or args.ids_file or "--ids"
    for name, user_labels, B in built:
        prefix = f"{args.prefix}_{name}" if name else args.prefix
        prefix = prefix.lstrip("_")
        title = args.title
        if name:
            title = f"{args.title} — {name}" if args.title else name
        n_rows = len(B["ids"])
        if n_rows == 0:
            print(f"{name or args.prefix}: nothing left to plot, skipped")
            continue
        print(f"{name or args.prefix}: {n_rows} rows "
              f"({B['stats'][0]} with a logFC or ON/OFF call in >=1 contrast, "
              f"{B['stats'][1]} quantified in >=1 pellet group)")
        print(f"  cells: {int(np.isfinite(B['diffmat']).sum())} logFC, "
              f"{int((B['onoff'] == 'ON').sum())} ON, {int((B['onoff'] == 'OFF').sum())} OFF, "
              f"{int(B['absent'].sum())} censored ({'/'.join(sorted(ABSENT_IN_MUTANTS))}), "
              f"{int((np.isnan(B['diffmat']) & (B['onoff'] == '') & ~B['absent']).sum())} ND")

        src_line = (f"# protein list: {src}"
                    + (f" (filters: {', '.join(args.where)})" if args.where else "")
                    + (f" [{args.split_sysid} system: {name}]" if name else "") + "\n")
        a = None
        if want_lfc:
            L = panel_view(B, args, "lfc")
            vlim = args.vlim or shared or finite_or(nanpct(L["diffmat"], 98), 1.0)
            plot_lfc(L["diffmat"], L["sig_cell"], L["onoff"], L["absent"], L["labels"],
                     L["blocks"], args, prefix, title, vlim, a)
            write_tsv(L, args, prefix, user_labels, order, "lfc", src_line, a)
        if want_abund:
            A = panel_view(B, args, "abundance")
            plot_abundance(A["mean"], A["in_matrix"], A["absent_group"], order,
                           A["labels"], A["blocks"], args, prefix, title, a)
            write_tsv(A, args, prefix, user_labels, order, "abundance", src_line, a)


def build_args():
    """Figure settings, fixed here rather than exposed as options."""
    sm = snakemake  # noqa: F821
    return SimpleNamespace(
        root=sm.params.root,
        mode="abundance",
        transpose=True,
        from_tsv=str(sm.input.candidates),
        prefix=sm.params.prefix,
        outdir=str(sm.params.outdir),
        gff=sm.params.gff,
        absent_replicon=sm.params.absent,
        abund_cmap="Spectral_r", abund_vmax=None, abund_vmin=None, annotate=False,
        cat_height=0.34, cmap="RdBu_r",
        data="DE/day2_log2_maxlfq_core_centred.tsv", day="2", drop_all_na=False,
        exclude=[], formats=["pdf", "svg"], group_col=None, id_col="MF6_ID",
        ids=None, ids_file=None, label_col="Description", label_fontsize=7.0,
        limma_dir="Limma", no_contig_last=False, no_filter=False, no_onoff=False,
        onoff_tier="weak", png=False, row_height=0.26, shared_vlim=False, sig=0.05,
        sig_col="adj_pval", sort="panel", split_sysid=None, title=None, vlim=None,
        where=[],
    )


if __name__ == "__main__":
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    main(build_args())
