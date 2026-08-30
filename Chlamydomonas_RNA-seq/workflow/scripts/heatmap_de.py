import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import matplotlib as mpl
mpl.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


ANNOTATION_COLS = [
    "Description", "Predalgo", "GO", "Best-hit-arabi-name"
]


def load_annotations(path):
    df = pd.read_csv(
        path, sep="\t", index_col="Gene",
        usecols=["Gene"] + ANNOTATION_COLS,
    )
    df.index.name = "gene"
    return df[ANNOTATION_COLS]


def load_clusters(up_path, down_path):
    usecols = ["Gene","Cluster"]
    up_df = pd.read_csv(up_path, usecols=usecols, sep="\t")
    up_df["Direction"] = "Up"
    down_df = pd.read_csv(down_path, usecols=usecols, sep="\t")
    down_df["Direction"] = "Down"

    df = pd.concat([up_df, down_df], ignore_index=True).set_index("Gene")

    df["Cluster_name"] = (
        df["Direction"] + "-" + df["Cluster"]
        .astype(int).astype(str)
    )
    df.index.name = "gene"
    df = (
        df.drop(columns=["Direction", "Cluster"])
        .rename(columns={"Cluster_name": "Cluster"})
    )

    return df


def load_data(deseq_path, subset_path):
    df_deseq = pd.read_csv(deseq_path, sep="\t")
    df_deseq = df_deseq.set_index("gene")

    df_subset = pd.read_csv(subset_path, sep="\t").dropna(subset=["gene"])
    df_subset = df_subset.drop_duplicates(subset="gene", keep="first")
    df_subset = df_subset.set_index("gene")

    df = df_subset.join(df_deseq)

    return df


def heatmap(df: pd.DataFrame,
            out_path,
            table_path=None,
            inch_x=8,
            inch_y=8,
            spacing=0.7,
            annotate=None,
            info_cols=None,
            col_as_name=None,
            sort_by=None,
            vmax=3,
            vmin=-3,
            info_cols_x=0.05,
            info_cols_y=-0.7,
            left_ticks=True):

    def padj_to_ast(x):
        if x <= 0.001:
            return "***"
        elif x <= 0.01:
            return "**"
        elif x <= 0.05:
            return "*"
        else:
            return ""
        
    if table_path is not None:
        # Build the table while the index is still the gene (Cre) ID, before
        # col_as_name overwrites it with the Alias below.
        lfc_tbl = df.pivot(columns="comparison", values="log2FoldChange")
        padj_tbl = df.pivot(columns="comparison", values="padj")
        stats_tbl = pd.DataFrame(index=lfc_tbl.index)
        for comp in lfc_tbl.columns:
            stats_tbl[f"{comp}_LFC"] = lfc_tbl[comp]
            stats_tbl[f"{comp}_adjPval"] = padj_tbl[comp]

        # Keep only Alias from the source gene list, plus the Mfuzz Cluster.
        keep_cols = [c for c in ["Alias", "Cluster"] if c in df.columns]
        info_tbl = df[keep_cols]
        info_tbl = info_tbl[~info_tbl.index.duplicated(keep="first")]

        table = info_tbl.join(stats_tbl).join(ANNOTATIONS)
        if sort_by in table.columns:
            table = table.sort_values(by=sort_by, kind="stable")
        table.index.name = "gene"
        table.to_csv(table_path, sep="\t")

    if col_as_name:
        new_index = df[col_as_name].fillna(pd.Series(df.index, index=df.index))
        df.index = new_index
    df.index.name = None

    df_pivot = df.pivot(columns="comparison",values="log2FoldChange").round(2)
    df_padj = df.pivot(columns="comparison",values="padj").map(padj_to_ast)
    df_annot = df_pivot.astype(str) + df_padj
    comparison_cols = df_pivot.columns
    df_pivot = df_pivot.join(
        df.drop(
            columns=["comparison", "log2FoldChange", "lfcSE", "pvalue", "padj"]
        )
    )
    if sort_by in df.columns:
        df_pivot = df_pivot.sort_values(by=sort_by, kind="stable")
    df_pivot = df_pivot.drop_duplicates(keep="first").fillna("")

    df_lfc = df_pivot[comparison_cols]
    df_lfc.columns = df_lfc.columns.str.split("_vs_").str[0]
    fig, ax = plt.subplots(figsize=(inch_x, inch_y))

    if annotate is True:
        annot = True
    elif annotate == "padj":
        annot = np.asarray(df_annot.loc[df_lfc.index])
    elif annotate == "*":
        annot = np.asarray(df_padj.loc[df_lfc.index])
    else:
        annot = False

    sns.heatmap(
        df_lfc, 
        ax=ax,
        cmap='RdBu_r',
        yticklabels=df_lfc.index,
        center=0,
        annot=annot,
        fmt="",
        vmax=vmax,
        vmin=vmin,
        linecolor="gray",
        cbar_kws={
            "location": "bottom", 
            "shrink": 1, "pad": 0.01, 
            "label": r"$\log_{2}\mathrm{FC}$"
        },
    )
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.tick_params(left=left_ticks, labelleft=left_ticks)
    tick_positions = ax.get_yticks()
    for c, info_col in enumerate(info_cols):
        ax.text(df_lfc.shape[1] + info_cols_x + (c * spacing), info_cols_y, info_col, va='center')
        for i, y in enumerate(tick_positions):
            ax.text(
                df_lfc.shape[1] + info_cols_x + (c * spacing), 
                y, 
                df_pivot.loc[df_lfc.index[i], info_col], 
                va='center'
            )

    plt.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close()


def make_plot(subset_file, output_file, table_file=None, **kwargs):
    df = load_data(deseq2_long, subset_file).join(UP_DOWN_CLUSTERS)
    df = df.dropna(subset="padj")
    heatmap(df, output_file, table_path=table_file, **kwargs)


if __name__ == "__main__":
    deseq2_long = snakemake.input.deseq2_long


    UP_DOWN_CLUSTERS = load_clusters(
        snakemake.input.up_clusters,
        snakemake.input.down_clusters
    )

    ANNOTATIONS = load_annotations(snakemake.input.all_annotations)
    
    make_plot(
        snakemake.input.gene_list,
        snakemake.output.heatmap,
        table_file=snakemake.output.get("table"),
        **snakemake.params.heatmap_arguments
    )


    