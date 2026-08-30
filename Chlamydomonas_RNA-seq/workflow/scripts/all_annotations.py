from snakemake.script import snakemake
import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd
import numpy as np

def merge_annotations(series):
    """
    Custom aggregator: removes NaNs, keeps unique values, 
    and concatenates them into a single string.
    """
    # Drop NaNs and convert to unique strings
    unique_values = series.dropna().astype(str).unique()

    if len(unique_values) == 0:
        return np.nan
    elif len(unique_values) == 1:
        return unique_values[0]
    else:
        # Split rows to words and join unique, keeping first-seen order.
        # A set here would reorder on every run (Python randomises string
        # hashing per process), so the output file - and every downstream
        # table built from it - differed between otherwise identical runs.
        all_words = dict.fromkeys(
            word for u in unique_values for word in u.split(" ")
        )
        return " ".join(all_words)

info_df = pd.read_csv(snakemake.input["annotation_info"], sep="\t")
# Remove unused cols
info_df = info_df.drop(columns=["#pacId", "transcriptName", "peptideName"])
# Group by locusName and apply the custom aggregator to all other columns
merged_df = info_df.groupby("locusName").agg(merge_annotations).reset_index()
# Remove _4532 from locus name
merged_df["locusName"] = merged_df["locusName"].str.split("_").str[0]
merged_df = merged_df.rename(columns={"locusName": "Gene"}).set_index("Gene")

master_df = pd.read_csv(snakemake.input["master_annotation"], sep="\t")
master_df["locusName_4532"] = master_df["locusName_4532"].str.split("_").str[0]
master_df = master_df.rename(columns={"locusName_4532": "Gene"}).set_index("Gene")
all_annotations = master_df.join(merged_df, how="outer")
# Del Removed genes e.g Cr_06_28453
all_annotations = all_annotations[all_annotations.index != ""]
all_annotations.to_csv(snakemake.output["all_annotations"], sep="\t")