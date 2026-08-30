from snakemake.script import snakemake
import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd

clusters_df = pd.read_csv(snakemake.input["clusters"], sep="\t", index_col=0)
membership_df = pd.read_csv(snakemake.input["membership"], sep="\t", index_col=0)
rename_clusters = snakemake.params.get("rename_clusters", None)
membership_df.columns = membership_df.columns.astype(int)

clusters_df['Membership_Score'] = [
    membership_df.loc[gene, cluster] 
    for gene, cluster in zip(clusters_df.index, clusters_df['Cluster'])
]

all_annotations_df = pd.read_csv(snakemake.input["all_annotations"], sep="\t", index_col=0)
clusters_info = clusters_df.join(all_annotations_df).sort_values(by="Membership_Score", kind="stable")
if rename_clusters is not None:
    rename_dict = rename_clusters.get(snakemake.wildcards.direction, {})
    clusters_info['Cluster'] = clusters_info['Cluster'].map(rename_dict)
clusters_info.to_csv(snakemake.output["clusters_annotated"], sep="\t")
