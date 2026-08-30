import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# # Load snakemake files and params
enrichment_table = snakemake.input["enrichment_table"]
enrichment_plot = snakemake.output["enrichment_plot"]

min_gene_count = snakemake.params["minGeneCount"]
min_pval = snakemake.params["minPval"]
# Plot all significant terms in QC plot

df = pd.read_csv(enrichment_table, sep="\t")
df_sig = df.loc[(df["Significant"] >= min_gene_count) 
                & (df["weight01Fisher"] < min_pval)]

# Filter top terms
df_top = df_sig
# Rich factor: the fraction of a term's annotated background genes that are
# in this gene set. NOT the conventional GeneRatio, which is
# Significant / (number of genes of interest).
df_top["RichFactor"] = df_top['Significant'] / df_top['Annotated']
df_top["Term_ont"] = df_top["Term"] + " (" + df_top["Ontology"] + ")"
fig, ax = plt.subplots(figsize=(10, len(df_top) * 0.25 + 2))

df_top = df_top.sort_values(by=["Ontology", "Cluster"])
# Scatterplot without legend
scatter = sns.scatterplot(
    data=df_top,
    x='RichFactor',
    y='Term_ont', 
    size='Significant',
    hue='neg_log10_p',
    sizes=(50, 400),
    palette='Spectral_r',
    edgecolor='black',
    hue_norm=(-np.log10(df_sig["weight01Fisher"].max()), 
              df_top['neg_log10_p'].max()), 
    linewidth=1,
    alpha=1,
    legend=False, 
    ax=ax
)

# Draw Legend colorbar
ax_cbar = ax.inset_axes([1.16, 0.50, 0.06, 0.4]) 

norm = plt.Normalize(-np.log10(df_sig["weight01Fisher"].max()), 
                     df_top['neg_log10_p'].max())
sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)
sm.set_array([])

cbar = fig.colorbar(sm, cax=ax_cbar, orientation='vertical')
ax_cbar.set_title(r"-log$_{10}$(P-value)", fontsize=10, pad=10) 

# Draw Legend for Gene Counts
ax_legend = ax.inset_axes([1.015, 0.45, 0, 0]) 
ax_legend.axis('off')

# Create manual dots
min_count = df_top['Significant'].min()
max_count = df_top['Significant'].max()
manual_sizes = [int(x) for x in np.linspace(min_count, max_count, 4)]

legend_elements = []
for size in manual_sizes:
    # Match the area scaling from scatterplot
    area = 50 + (size - min_count) / (max_count - min_count) * (400 - 50)
    marker_diam = np.sqrt(area) 
    
    legend_elements.append(
        Line2D([0], [0], marker='o', color='w', label=str(size),
               markerfacecolor='#D6D6D6', markersize=marker_diam, 
               markeredgecolor='black')
    )

# Add legend
ax_legend.legend(handles=legend_elements, loc='upper left', 
                 title="Gene Count", frameon=False, labelspacing=0.8)

ax.grid(True, linestyle='--', alpha=0.3)
ax.set_xlabel("Rich factor", fontsize=12)
ax.set_ylabel("GO Term", fontsize=12)
plt.tight_layout()
plt.savefig(enrichment_plot, dpi=300, bbox_inches='tight')