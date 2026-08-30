import sys
from pathlib import Path

import pandas as pd
from snakemake.script import snakemake

# send stdout and stderr to the rule's log file
sys.stdout = open(snakemake.log[0], "w")
sys.stderr = sys.stdout

column_name = snakemake.params["column_name"]

dfs = []
for file in snakemake.input:
    file = Path(file)
    df = pd.read_csv(file, sep='\t', index_col="gene_id")[[column_name]]
    df.columns = [file.stem.split('.')[0]]
    dfs.append(df)

matrix = pd.concat(dfs, axis=1)
matrix.index.name = "gene"
matrix.index = matrix.index.str.split("_").str[0]
matrix = matrix.sort_index(axis=1).sort_index()
matrix.to_csv(snakemake.output[0], sep='\t')
