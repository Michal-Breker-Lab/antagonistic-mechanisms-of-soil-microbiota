#!/usr/bin/env python3
"""PCA of the log2 MaxLFQ matrix, one figure per sample selection."""

import sys
from types import SimpleNamespace
import os
import re

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
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
import seaborn as sns

STRAIN_ORDER = ["MF6", "27D6", "34F7"]
DAY_MARK = {"2": "^", "3": "o"}
META_COLS = ["Protein", "Gene", "Protein Length", "Combined Total Peptides",
             "Combined Unique Spectral Count", "Protein Probability", "Description"]
SAMPLE_RE = re.compile(r"^(?P<strain>[^_]+)_(?P<cond>alone|withC)_d(?P<day>[23])_"
                       r"(?P<rep>[1-4])$")
QC_FAILED = []
ABSENT_IN_MUTANTS = set()


def parse_args():
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")  # noqa: F821
    sm = snakemake  # noqa: F821
    args = SimpleNamespace(
        root=sm.params.root,
        data=sm.params.data,
        gff=sm.params.gff,
        outdir=sm.params.outdir,
        prefix=sm.wildcards.variant,
        days=list(sm.params.days),
        cond=list(sm.params.cond),
        exclude=list(sm.params.exclude),
        absent_replicon=sm.params.absent,
        pcs=6,
        png=False,
    )
    global ABSENT_IN_MUTANTS, QC_FAILED
    ABSENT_IN_MUTANTS = {args.absent_replicon}
    QC_FAILED = list(sm.params.qc_failed)
    return args


def contig_of(path):
    out = {}
    if not os.path.exists(path):
        return out
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


def eta_sq(scores, groups):
    """Share of a PC's variance that lies BETWEEN the levels of one factor.

    1.0 = the factor explains the axis completely, 0 = not at all.  This is what
    turns "PC1 is 34% of the variance" into "PC1 is the day axis"."""
    scores, groups = np.asarray(scores, float), np.asarray(groups)
    grand = scores.mean()
    ss_tot = ((scores - grand) ** 2).sum()
    if ss_tot == 0:
        return np.nan
    ss_between = sum(len(scores[groups == g]) * (scores[groups == g].mean() - grand) ** 2
                     for g in np.unique(groups))
    return ss_between / ss_tot


def main():
    args = parse_args()
    outdir = os.path.join(args.root, args.outdir)
    os.makedirs(outdir, exist_ok=True)

    df = pd.read_csv(os.path.join(args.root, args.data), sep="\t")
    drop, days, conds = set(args.exclude), set(args.days), set(args.cond)
    samples = [c for c in df.columns
               if c not in META_COLS and SAMPLE_RE.match(c) and c not in drop
               and SAMPLE_RE.match(c).group("day") in days
               and SAMPLE_RE.match(c).group("cond") in conds]
    if not samples:
        raise SystemExit(f"no runs left for day(s) {'/'.join(sorted(days))} and "
                         f"condition(s) {'/'.join(sorted(conds))} after exclusions")
    meta = pd.DataFrame([SAMPLE_RE.match(s).groupdict() for s in samples], index=samples)

    X = df[samples].apply(pd.to_numeric, errors="coerce").replace(0, np.nan)
    contig = contig_of(os.path.join(args.root, args.gff))
    on_c2 = df["Protein"].map(contig).isin(ABSENT_IN_MUTANTS)
    complete = (~X.isna()).all(axis=1)
    core = complete & ~on_c2
    n_art = int((complete & on_c2).sum())
    L = np.log2(X[core].to_numpy(float))
    print(f"{len(samples)} runs x {len(df)} proteins -> core set of {int(core.sum())} "
          f"quantified in every run"
          + (f"; {n_art} further complete row(s) censored as "
             f"{'/'.join(sorted(ABSENT_IN_MUTANTS))} match-between-runs artifacts"
             if n_art else ""))
    applied = sorted(c for c in drop
                     if (m := SAMPLE_RE.match(c)) and m.group("day") in days
                     and m.group("cond") in conds)
    if applied:
        print(f"excluded: {', '.join(applied)}")

    L = L - np.median(L, axis=0) + np.median(np.median(L, axis=0))
    M = (L - L.mean(axis=1, keepdims=True)).T

    U, S, _ = np.linalg.svd(M, full_matrices=False)
    scores = U * S
    var = S ** 2 / (S ** 2).sum() * 100
    npc = min(args.pcs, scores.shape[1])

    fac = {k: v for k, v in (("strain", meta["strain"].to_numpy()),
                             ("Chlamy", meta["cond"].to_numpy()),
                             ("day", meta["day"].to_numpy()))
           if len(set(v)) > 1}
    print("\nvariance explained, and what each axis separates (eta^2, 1 = the axis IS "
          "that factor):")
    print(f"  {'PC':>4} {'var%':>7}  " + "  ".join(f"{k:>7}" for k in fac))
    for i in range(npc):
        e = {k: eta_sq(scores[:, i], v) for k, v in fac.items()}
        print(f"  {'PC'+str(i+1):>4} {var[i]:7.1f}  "
              + "  ".join(f"{e[k]:7.2f}" for k in fac))

    sc = pd.DataFrame(scores[:, :npc], index=samples,
                      columns=[f"PC{i+1}" for i in range(npc)])
    sc.insert(0, "qc_failed", ["yes" if s in QC_FAILED else "no" for s in samples])
    for c in ("strain", "cond", "day", "rep"):
        sc.insert(1, c, meta[c])
    sc.index.name = "run"
    with open(os.path.join(outdir, f"{args.prefix}_scores.tsv"), "w") as fh:
        sel = (f"day {' + '.join(sorted(days))}; "
               f"{' + '.join(sorted(conds))} culture")
        fh.write(f"# PCA of {len(samples)} runs over {int(core.sum())} proteins "
                 f"quantified in every one of them ({sel})\n"
                 f"# source: {args.data}, log2, core-set median centred per run\n"
                 "# qc_failed = yes -> shallow AND contaminant-dominated in the QC report "
                 "(only present if --exclude was overridden)\n")
        sc.to_csv(fh, sep="\t")

    va = pd.DataFrame({"PC": [f"PC{i+1}" for i in range(npc)],
                       "variance_pct": var[:npc].round(3)})
    for k, v in fac.items():
        va[f"eta2_{k}"] = [round(eta_sq(scores[:, i], v), 4) for i in range(npc)]
    with open(os.path.join(outdir, f"{args.prefix}_variance.tsv"), "w") as fh:
        fh.write("# variance explained per PC, and how much of each PC's variance lies "
                 "between the levels of a factor\n"
                 "# eta2 = between-group SS / total SS: 1.0 = that PC is exactly that "
                 "factor, 0 = unrelated to it\n")
        va.to_csv(fh, sep="\t", index=False)

    pal = dict(zip(STRAIN_ORDER, sns.color_palette("Set2", 3)))
    fig = plt.figure(figsize=(9.6, 6.2))
    gs = GridSpec(2, 2, width_ratios=[3.05, 1.0], height_ratios=[1.25, 1.0],
                  wspace=0.06, hspace=0.55, figure=fig)
    ax = fig.add_subplot(gs[:, 0])
    axl = fig.add_subplot(gs[0, 1]); axl.axis("off")
    axs = fig.add_subplot(gs[1, 1])

    for i, s in enumerate(samples):
        m = meta.loc[s]
        filled = m["cond"] == "withC"
        col = pal.get(m["strain"], "#777777")
        flagged = s in QC_FAILED
        ax.plot(scores[i, 0], scores[i, 1], marker=DAY_MARK.get(m["day"], "o"),
                markersize=8.5, linestyle="none",
                markerfacecolor=col if filled else "none", markeredgecolor=col,
                markeredgewidth=1.4, zorder=3)
        if flagged:
            ax.plot(scores[i, 0], scores[i, 1], marker="o", markersize=15.0,
                    linestyle="none", markerfacecolor="none", markeredgecolor="black",
                    markeredgewidth=1.1, zorder=4)
            ax.annotate(f"{s}  (QC-failed)", (scores[i, 0], scores[i, 1]),
                        textcoords="offset points", xytext=(12, -3), fontsize=7,
                        color="#333333")
    ax.axhline(0, color="#cccccc", lw=0.8, zorder=1)
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=1)
    ax.set_xlabel(f"PC1  ({var[0]:.1f}% of variance)", fontsize=10)
    ax.set_ylabel(f"PC2  ({var[1]:.1f}%)", fontsize=10)
    ax.tick_params(labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("black"); sp.set_linewidth(0.8)

    keys = [Line2D([], [], marker="o", linestyle="none", markersize=8.5,
                   markerfacecolor=pal[s], markeredgecolor=pal[s], label=s)
            for s in STRAIN_ORDER if s in set(meta["strain"])]
    present = set(meta["cond"])
    keys += [Line2D([], [], marker="o", linestyle="none", markersize=8.5,
                    markerfacecolor="#555555" if c == "withC" else "none",
                    markeredgecolor="#555555",
                    label=("+ Chlamydomonas" if c == "withC" else "− Chlamydomonas"))
             for c in ("withC", "alone") if c in present]
    keys += [Line2D([], [], marker=DAY_MARK[d], linestyle="none", markersize=8.5,
                    markerfacecolor="none", markeredgecolor="#555555", label=f"day {d}")
             for d in sorted(set(meta["day"]))]
    axl.legend(handles=keys, loc="upper left", bbox_to_anchor=(-0.06, 1.02),
               frameon=False, fontsize=8.5, handletextpad=0.7, labelspacing=0.75)

    axs.bar(range(1, npc + 1), var[:npc], color="#b9b8b2", edgecolor="black", lw=0.6)
    axs.set_xticks(range(1, npc + 1))
    axs.set_xticklabels([str(i) for i in range(1, npc + 1)], fontsize=7.5)
    axs.set_xlabel("principal component", fontsize=8)
    axs.set_ylabel("% variance", fontsize=8)
    axs.tick_params(labelsize=7.5)
    for sp in ("top", "right"):
        axs.spines[sp].set_visible(False)

    for ext in ("pdf", "svg"):
        out = os.path.join(outdir, f"{args.prefix}.{ext}")
        fig.savefig(out, bbox_inches="tight")
        print(f"wrote {out}")
    if args.png:
        out = os.path.join(outdir, f"{args.prefix}.png")
        fig.savefig(out, bbox_inches="tight", dpi=150)
        print(f"wrote {out}")
    plt.close(fig)
    print(f"wrote {outdir}/{args.prefix}_scores.tsv\n"
          f"wrote {outdir}/{args.prefix}_variance.tsv")


if __name__ == "__main__":
    main()
