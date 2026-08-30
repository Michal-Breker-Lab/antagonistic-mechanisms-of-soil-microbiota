from snakemake.script import snakemake
import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd

use_clusters = snakemake.params["use_clusters"]

cluster_df = pd.read_csv(snakemake.input["clusters"], sep="\t")
cluster_df = cluster_df.loc[cluster_df["Cluster"].isin(use_clusters)]
cluster_df["Cluster"] = "-".join(list(map(str,use_clusters)))
cluster_df["Membership_Score"] = 1

cluster_df.to_csv(snakemake.output["clusters"], index=False, sep="\t")
