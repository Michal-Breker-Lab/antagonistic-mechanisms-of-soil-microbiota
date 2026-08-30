"""Horizontal (transposed) version of heatmap_de.py.

Genes run along the x-axis and the Ts-vs-T0 comparisons along the y-axis; the
Mfuzz cluster is a bracketed group label under the heatmap instead of a text
column on the right.

The figure is written as PDF with real, editable text (pdf.fonttype 42 embeds a
subsetted TrueType font instead of matplotlib's default Type-3), so Illustrator
and Acrobat see selectable text rather than outlines.

Layout is driven by params.heatmap_arguments; see the heatmap() signature.
"""

import sys

sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")

import matplotlib as mpl

mpl.use("Agg")

# Keep text as text in the PDF, in a single sans family.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
mpl.rcParams["mathtext.fontset"] = "custom"
for _key in ("mathtext.rm", "mathtext.it", "mathtext.bf", "mathtext.sf"):
    mpl.rcParams[_key] = mpl.rcParams["font.sans-serif"][0]
mpl.rcParams["mathtext.default"] = "regular"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.transforms import offset_copy


def padj_to_ast(x):
    if x <= 0.001:
        return "***"
    if x <= 0.01:
        return "**"
    if x <= 0.05:
        return "*"
    return ""


def load_clusters(up_path, down_path):
    usecols = ["Gene", "Cluster"]
    up_df = pd.read_csv(up_path, usecols=usecols, sep="\t")
    up_df["Direction"] = "Up"
    down_df = pd.read_csv(down_path, usecols=usecols, sep="\t")
    down_df["Direction"] = "Down"

    df = pd.concat([up_df, down_df], ignore_index=True).set_index("Gene")
    df["Cluster_name"] = df["Direction"] + "-" + df["Cluster"].astype(int).astype(str)
    df.index.name = "gene"
    df = df.drop(columns=["Direction", "Cluster"]).rename(
        columns={"Cluster_name": "Cluster"}
    )
    return df


def load_data(deseq_path, subset_path):
    df_deseq = pd.read_csv(deseq_path, sep="\t").set_index("gene")

    df_subset = pd.read_csv(subset_path, sep="\t").dropna(subset=["gene"])
    df_subset = df_subset.drop_duplicates(subset="gene", keep="first")
    df_subset = df_subset.set_index("gene")

    return df_subset.join(df_deseq)


def prepare(df, col_as_name=None, sort_by=None, annotate="padj"):
    """Reshape exactly as heatmap_de.py does, so gene order and values match.

    Returns (lfc, labels, stars, cluster) with genes still on the index; the
    caller transposes. Values and significance stars are kept apart so they can
    be placed independently inside a cell.
    """
    if col_as_name:
        df = df.copy()
        df.index = df[col_as_name].fillna(pd.Series(df.index, index=df.index))
    df.index.name = None

    df_pivot = df.pivot(columns="comparison", values="log2FoldChange").round(2)
    df_padj = df.pivot(columns="comparison", values="padj").map(padj_to_ast)
    df_label = df_pivot.astype(str)

    comparison_cols = df_pivot.columns
    df_pivot = df_pivot.join(
        df.drop(columns=["comparison", "log2FoldChange", "lfcSE", "pvalue", "padj"])
    )
    if sort_by in df.columns:
        df_pivot = df_pivot.sort_values(by=sort_by, kind="stable")
    df_pivot = df_pivot.drop_duplicates(keep="first").fillna("")

    lfc = df_pivot[comparison_cols]
    lfc.columns = lfc.columns.str.split("_vs_").str[0]

    labels = df_label.loc[lfc.index] if annotate in ("padj", True) else None
    stars = df_padj.loc[lfc.index] if annotate in ("padj", "*") else None

    cluster = df_pivot["Cluster"] if "Cluster" in df_pivot.columns else None
    return lfc, labels, stars, cluster


def relative_luminance(colors):
    """sRGB relative luminance - the rule seaborn uses to pick text colour."""
    rgb = np.asarray(colors)[:, :3]
    rgb = np.where(rgb <= 0.03928, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    return rgb @ np.array([0.2126, 0.7152, 0.0722])


def text_width(ax, s, fontsize, renderer):
    """Width of `s` in display pixels at the annotation font size."""
    probe = ax.text(0, 0, s, fontsize=fontsize)
    width = probe.get_window_extent(renderer=renderer).width
    probe.remove()
    return width


def annotate_cells(fig, ax, labels, stars, annot_size, num_dy, star_dy):
    """Draw log2FC values and significance stars inside the cells.

    Done by hand rather than with seaborn's `annot=` so that the digits are
    centred on the cell: a centred "-0.42" would put its *string* centre on the
    cell centre, leaving the digits half a minus sign too far right. Negative
    values are therefore nudged left by half a minus sign, applied in points at
    draw time so it stays exact after the equal-aspect rescaling of the axes.
    """
    mesh = ax.collections[0]
    mesh.update_scalarmappable()
    masked = np.ma.getmaskarray(mesh.get_array()).ravel()
    lum = relative_luminance(mesh.get_facecolors())

    renderer = fig.canvas.get_renderer()
    minus_px = text_width(ax, "0-", annot_size, renderer) - text_width(
        ax, "0", annot_size, renderer
    )
    neg_transform = offset_copy(
        ax.transData, fig=fig, x=-minus_px / 2 / fig.dpi * 72, y=0, units="points"
    )

    n_rows, n_cols = (labels if labels is not None else stars).shape
    for k, (i, j) in enumerate(np.ndindex(n_rows, n_cols)):
        if masked.size and masked[k]:
            continue
        color = ".15" if lum[k] > 0.408 else "w"
        x, y = j + 0.5, i + 0.5

        if labels is not None:
            label = str(labels.iat[i, j])
            ax.text(
                x,
                y + num_dy,
                label,
                transform=neg_transform if label.startswith("-") else ax.transData,
                ha="center",
                va="center",
                fontsize=annot_size,
                color=color,
            )
        if stars is not None and str(stars.iat[i, j]):
            ax.text(
                x,
                y + (star_dy if labels is not None else 0),
                str(stars.iat[i, j]),
                ha="center",
                va="center",
                fontsize=annot_size,
                color=color,
            )


def cluster_blocks(cluster):
    """Contiguous runs of the (already sorted) cluster column -> brackets."""
    blocks = []
    start = 0
    values = list(cluster.fillna("").astype(str))
    for i in range(1, len(values) + 1):
        if i == len(values) or values[i] != values[start]:
            blocks.append((start, i, values[start]))
            start = i
    return blocks


def heatmap(
    df,
    out_path,
    col_as_name=None,
    sort_by=None,
    annotate="padj",
    vmin=-2,
    vmax=2,
    cell_w=0.52,
    cell_h=0.52,
    annot_size=7,
    tick_size=8,
    cluster_size=8,
    gene_rotation=0,
    num_dy=0.0,
    star_dy=0.22,
    bracket_y=0.22,
    bracket_label_y=0.62,
    cluster_label="Cluster",
    cbar_label=r"$\log_{2}$FC",
):
    lfc, labels, stars, cluster = prepare(
        df, col_as_name=col_as_name, sort_by=sort_by, annotate=annotate
    )

    # Transpose: genes -> columns, T1..T4 -> rows.
    mat = lfc.T
    labels = labels.T if labels is not None else None
    stars = stars.T if stars is not None else None

    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(
        figsize=(n_cols * cell_w + 2.4, n_rows * cell_h + 2.6)
    )

    sns.heatmap(
        mat,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=vmin,
        vmax=vmax,
        linewidths=0,  # solid heatmap: no gaps between cells
        square=True,
        cbar_kws={
            "location": "right",
            "shrink": 0.55,
            "pad": 0.02,
            "aspect": 10,
            "label": cbar_label,
        },
    )
    ax.set_xlabel("")
    ax.set_ylabel("")

    if labels is not None or stars is not None:
        annotate_cells(fig, ax, labels, stars, annot_size, num_dy, star_dy)

    # Gene names on top, comparisons on the left.
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    plt.setp(
        ax.get_xticklabels(),
        rotation=gene_rotation,
        ha="center" if gene_rotation in (0, 90) else "left",
        va="bottom",
        fontsize=tick_size,
    )
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=tick_size)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=tick_size)
    cbar.set_label(cbar_label, fontsize=tick_size)

    # Mfuzz clusters as brackets under the heatmap. Genes are sorted by
    # cluster, so every cluster is one contiguous block of columns.
    if cluster is not None:
        y0 = n_rows + bracket_y
        y1 = n_rows + bracket_label_y
        for start, stop, name in cluster_blocks(cluster):
            if not name:
                continue
            ax.plot(
                [start + 0.12, stop - 0.12],
                [y0, y0],
                color="black",
                lw=0.8,
                clip_on=False,
                solid_capstyle="butt",
            )
            ax.text(
                (start + stop) / 2,
                y1,
                name,
                ha="center",
                va="center",
                fontsize=cluster_size,
                clip_on=False,
            )
        ax.text(
            -0.25,
            n_rows + (bracket_y + bracket_label_y) / 2,
            cluster_label,
            ha="right",
            va="center",
            fontsize=cluster_size,
            clip_on=False,
        )

    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    up_down_clusters = load_clusters(
        snakemake.input.up_clusters, snakemake.input.down_clusters
    )

    data = load_data(snakemake.input.deseq2_long, snakemake.input.gene_list).join(
        up_down_clusters
    )
    data = data.dropna(subset="padj")

    heatmap(data, snakemake.output.heatmap, **snakemake.params.heatmap_arguments)
