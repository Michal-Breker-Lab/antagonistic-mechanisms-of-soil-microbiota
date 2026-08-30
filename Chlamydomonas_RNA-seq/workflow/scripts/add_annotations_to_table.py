from snakemake.script import snakemake
import sys
sys.stderr = sys.stdout
sys.stdout = open(snakemake.log[0], "w")
import pandas as pd

ANNOTATION_COLS = [
    "geneSymbol", "Description", "Predalgo", "GO", "Best-hit-arabi-name"
]

table = pd.read_csv(snakemake.input["table"], sep="\t", index_col="gene")

annotations = pd.read_csv(
    snakemake.input["all_annotations"], sep="\t", index_col="Gene",
    usecols=["Gene"] + ANNOTATION_COLS,
)
annotations.index.name = "gene"

table = table.join(annotations[ANNOTATION_COLS])
table.to_csv(snakemake.output["annotated"], sep="\t")
