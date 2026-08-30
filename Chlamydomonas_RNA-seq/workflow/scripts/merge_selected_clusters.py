from snakemake.script import snakemake
import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd

clusters_up = pd.read_csv(snakemake.input["up_clusters"], sep="\t")
clusters_down = pd.read_csv(snakemake.input["down_clusters"], sep="\t")

selected_up = snakemake.params["up_clusters"]
selected_down = snakemake.params["down_clusters"]
clusters_up_selected = clusters_up.loc[clusters_up["Cluster"].isin(selected_up)]
clusters_down_selected = clusters_down.loc[clusters_down["Cluster"].isin(selected_down)]

merged_clusters = pd.concat([clusters_up_selected, clusters_down_selected], ignore_index=True)
merged_clusters["Cluster"] = "UP-" + "-".join(list(map(str,selected_up))) + "_DOWN-" + "-".join(list(map(str,selected_down)))
merged_clusters["Membership_Score"] = 1
merged_clusters.to_csv(snakemake.output["merged_table"], sep="\t", index=False)