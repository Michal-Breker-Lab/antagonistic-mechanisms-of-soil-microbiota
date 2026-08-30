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

zscores = pd.read_csv(snakemake.input["zscores"], sep="\t", index_col=0)
clusters_df = pd.read_csv(snakemake.input["cluster_table"], sep="\t", index_col=0)

# Filter low membership genes in QC plots
min_membership = float(snakemake.params.get("minMembership", 0))
if min_membership is not None:
    clusters_df = clusters_df[clusters_df["Membership_Score"] >= min_membership]

# Sort zscores based on cluster assignment so the heatmap rows match the plots
zscores_sorted = zscores.loc[clusters_df.sort_values(by="Cluster", kind="stable").index]

# 2. Setup Figure Layout
n_clusters = clusters_df["Cluster"].max()
fig = plt.figure(figsize=(10, 2 * n_clusters))
gs = GridSpec(n_clusters, 2, width_ratios=[0.3, 1], wspace=0.01, hspace=0.05)

cmap_mems = plt.get_cmap('Spectral_r')

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
    
    # Clean up heatmap x-axis
    if i < n_clusters - 1:
        ax_heat.set_xticks([])
        ax_line.set_xticks([])

# Add legend to Line Plot
# Get coordinates of the first line plot axis
pos_line = fig.axes[0].get_position() 
left_line = pos_line.x0
width_line = pos_line.x1 - pos_line.x0
legend_width = 0.15
# Center the legend over this width
cax1 = fig.add_axes([left_line + (width_line - legend_width)/2, 0.9, legend_width, 0.02]) 

norm1 = mpl.colors.Normalize(vmin=0, vmax=1)
cb1 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm1, cmap=cmap_mems),
                   cax=cax1, orientation='horizontal', ticks=[0.2, 0.4, 0.6, 0.8])
cb1.ax.set_title('Membership', pad=5, fontsize=10)

# Add legend to Heatmap
pos_heat = fig.axes[1].get_position()
left_heat = pos_heat.x0
width_heat = pos_heat.x1 - pos_heat.x0
# Center the legend over the heatmap column
cax2 = fig.add_axes([left_heat + (width_heat - legend_width)/2, 0.9, legend_width, 0.02])
norm2 = mpl.colors.Normalize(vmin=-2, vmax=2)
cb2 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm2, cmap="RdBu_r"),
                   cax=cax2, orientation='horizontal', ticks=[-2, -1, 0, 1, 2])
cb2.ax.set_title('Expression', pad=5, fontsize=10)

# Final spacing adjustment
plt.subplots_adjust(top=0.88)
plt.savefig(snakemake.output["expression_plot"], bbox_inches="tight")

