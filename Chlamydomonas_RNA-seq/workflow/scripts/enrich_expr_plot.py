from snakemake.script import snakemake
import sys

sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")

import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import numpy as np
import matplotlib.lines as mlines

def split_long_term(text):
    # Check if text is longer than 20 chars
    if len(str(text)) > 20:
        words = str(text).split()
        # Ensure there are at least 3 words to split between
        if len(words) >= 3:
            # Reconstruct: "Word1 Word2\nWord3..."
            return f"{words[0]} {words[1]}\n{' '.join(words[2:])}"
        elif len(words[0]) > 15:
            return f"{words[0]}\n{' '.join(words[1:])}"
    return text

 
zscores = pd.read_csv(snakemake.input["zscores"], sep="\t", index_col=0)
clusters_df = pd.read_csv(snakemake.input["cluster_table"], sep="\t", index_col=0)
enrichment_table = snakemake.input["enrichment_table"]
min_gene_count = snakemake.params["minGeneCount"]
min_pval = snakemake.params["minPval"]
# Plot all significant terms in QC plot
top_n = snakemake.params.get("topNterms", 1000)

df_enrich = pd.read_csv(enrichment_table, sep="\t")
df_sig = df_enrich.loc[(df_enrich["Significant"] >= min_gene_count) 
                & (df_enrich["weight01Fisher"] < min_pval)]
# Filter top terms
df_top = df_sig.groupby(['Cluster']).apply(
    lambda x: x.nlargest(top_n, 'neg_log10_p')
).reset_index(drop=False)
df_top['Term'] = df_top['Term'].apply(split_long_term)
df_top = df_top.sort_values(by="neg_log10_p", ascending=True, kind="stable")

ontology_order = sorted(df_top["Ontology"].unique())
ontology_map = {val: i for i, val in enumerate(ontology_order)}
min_count = df_top['Significant'].min()
max_count = df_top['Significant'].max()
min_logpval = -np.log10(min_pval)
max_logpval = round(df_top['neg_log10_p'].max())
size_range = (50, 400)

# Filter low membership genes in QC plots
min_membership = float(snakemake.params.get("minMembership", 0))
if min_membership is not None:
    clusters_df = clusters_df[clusters_df["Membership_Score"] >= min_membership]

# 2. Setup Figure Layout
n_clusters = clusters_df["Cluster"].max()
fig = plt.figure(figsize=(12, 2.5 * n_clusters))
gs = GridSpec(n_clusters, 3, width_ratios=[0.3, 1, 0.3], wspace=0.01, hspace=0.05)

cmap_mems = plt.get_cmap('Spectral_r')
cmap_enrich = plt.get_cmap('Spectral_r')
for i in range(n_clusters):
    cluster_id = i + 1
    
    # Mfuzz Trend Plot
    ax_line = fig.add_subplot(gs[i, 0])
    cluster_df = clusters_df[clusters_df["Cluster"] == cluster_id]
    cluster_genes = cluster_df.index
    mems = cluster_df["Membership_Score"].sort_values(ascending=True, kind="stable")
    # Plot gene lines
    for gene in mems.index:
        ax_line.plot(zscores.columns, zscores.loc[gene], 
                     color=cmap_mems(mems.loc[gene]), alpha=1,
                     lw=0.5)

    ax_line.set_ylabel(f"Cluster {cluster_id}\nn={len(cluster_genes)}", 
                       rotation=0, va='center',ha='right', fontweight='bold',
                       multialignment='center')
    ax_line.set_yticks([])
    for _, spine in ax_line.spines.items():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1)

    # Heatmap Segment
    ax_heat = fig.add_subplot(gs[i, 1])
    sns.heatmap(zscores.loc[cluster_genes], 
                ax=ax_heat, 
                cmap="RdBu_r", 
                center=0, 
                vmin=-2,
                vmax=2,
                cbar=False,
                yticklabels=False)
    ax_heat.set_ylabel("")
    for _, spine in ax_heat.spines.items():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1)
    
    ax_enrich = fig.add_subplot(gs[i, 2])
    enrich_cn = f"{snakemake.wildcards['direction'].capitalize()}_C{cluster_id}"
    df_top_sub = df_top[df_top["Cluster"] == enrich_cn].copy()
    df_top_sub['Ontology_Idx'] = df_top_sub['Ontology'].map(ontology_map)
    scatter = sns.scatterplot(
        data=df_top_sub,
        x='Ontology_Idx',
        y='Term',
        size='Significant',
        sizes=size_range,
        size_norm=(min_count, max_count),
        hue="neg_log10_p",
        palette=cmap_enrich,
        hue_norm=(min_logpval, max_logpval), 
        edgecolor='black',
        linewidth=1,
        alpha=1,
        legend=False, 
        ax=ax_enrich
    )

    ax_enrich.set_title(None)

    ax_enrich.set_xticks(range(len(ontology_order)))
    ax_enrich.set_xticklabels(ontology_order)
    ax_enrich.set_xlim(-0.5, 2.5)
    ax_enrich.set_xlabel(None)

    ax_enrich.set_ylabel(None)
    ax_enrich.yaxis.set_ticks_position('right')
    ax_enrich.yaxis.set_label_position('right')
    ax_enrich.set_ylim(-0.5, max(df_top_sub.shape[0] - 0.5, 0.5)) 

    ax_enrich.grid(True, linestyle='--', alpha=0.3)

    if df_top_sub.empty:
        ax_enrich.tick_params(axis='y', right=False, labelright=False) 

    # Clean up x-axis
    if i < n_clusters - 1:
        ax_heat.set_xticks([])
        ax_line.set_xticks([])
        ax_enrich.tick_params(axis='x', bottom=False) 
    
pad = 8
legend_width = 0.15
height = 0.02
y_pos = 0.9
# Add legend to Line Plot
# Get coordinates of the first line plot axis
pos_line = fig.axes[0].get_position() 
left_line = pos_line.x0
width_line = pos_line.x1 - pos_line.x0
# Center the legend over this width
cax1 = fig.add_axes([left_line + (width_line - legend_width)/2, y_pos, legend_width, height]) 
norm1 = mpl.colors.Normalize(vmin=0, vmax=1)
cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm1, cmap=cmap_mems),
                   cax=cax1, orientation='horizontal', ticks=[0.2, 0.4, 0.6, 0.8])
cb1.ax.set_title('Membership', pad=pad, fontsize=10)

# Add legend to Heatmap
pos_heat = fig.axes[1].get_position()
left_heat = pos_heat.x0
width_heat = pos_heat.x1 - pos_heat.x0
# Center the legend over the heatmap column
cax2 = fig.add_axes([left_heat + (width_heat - legend_width)/2, y_pos, legend_width, height])
norm2 = mpl.colors.Normalize(vmin=-2, vmax=2)
cb2 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm2, cmap="RdBu_r"),
                   cax=cax2, orientation='horizontal', ticks=[-2, -1, 0, 1, 2])
cb2.ax.set_title('Expression', pad=pad, fontsize=10)

# Add legend to Enrichment Plot
pos_enrich = fig.axes[2].get_position() 
left_enrich = pos_enrich.x0
width_enrich = pos_enrich.x1 - pos_enrich.x0
# Calculate starting X to center the pair (colorbar + gap + size_legend)
start_x = left_enrich + (width_enrich - legend_width) / 2
# Colorbar Axis (Left)
cax1 = fig.add_axes([start_x, y_pos, legend_width, height]) 
norm1 = mpl.colors.Normalize(vmin=min_logpval, vmax=max_logpval)
ticks = list(np.linspace(min_logpval, max_logpval, 5, dtype=int))
cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm1, cmap=cmap_enrich),
                   cax=cax1, orientation='horizontal', ticks=ticks)
cb1.ax.set_title(r"-log$_{10}$(P-value)", pad=pad, fontsize=10)

# Create a dedicated axis for the dots
cax2 = fig.add_axes([start_x + legend_width + 0.01, y_pos, legend_width, height])
cax2.set_axis_off()  # Hide borders/ticks
cax2.set_title("Gene Count", pad=pad, fontsize=10) # Title on the axis itself
# Define sizes
manual_sizes = [int(x) for x in np.linspace(min_count, max_count, 4)]
# Create the legend on this new axis
legend_elements = []
for size in manual_sizes:
    # Match the area scaling from scatterplot
    area = (size_range[0] + (size - min_count) 
            / (max_count - min_count) 
            * (size_range[1] - size_range[0]))
    marker_diam = np.sqrt(area)
    
    legend_elements.append(
        mlines.Line2D([], [], marker='o', color='w', label=str(size),
                      markerfacecolor='#D6D6D6', markersize=marker_diam, 
                      markeredgecolor='black')
    )

# Add the legend to the custom axis
reordered_elements = [legend_elements[0], legend_elements[2], 
                      legend_elements[1], legend_elements[3]]
cax2.legend(handles=reordered_elements, loc='center', 
            ncol=2, frameon=False, 
            handletextpad=0.2, 
            columnspacing=0.5,
            labelspacing=0.8)

# Final spacing adjustment
plt.subplots_adjust(top=0.88)
plt.savefig(snakemake.output["expression_plot"], bbox_inches="tight")

