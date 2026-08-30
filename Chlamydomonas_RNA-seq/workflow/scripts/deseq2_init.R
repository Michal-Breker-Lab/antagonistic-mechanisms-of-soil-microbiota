log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)

# Count Data
counts_data <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
counts_data <- counts_data[, order(colnames(counts_data))]

keep <- rowSums(counts_data >= 10) >= 6
counts_filtered <- counts_data[keep, ]

# Sample Metadata
col_data <- read.table(
  snakemake@config[["samples"]],
  sep = ",",
  header = TRUE,
  row.names = "sample_name",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

condition_col <- snakemake@params[["condition"]]
col_data[[condition_col]] <- as.factor(col_data[[condition_col]])

# DESeq2 Object
design_formula <- as.formula(paste0("~", condition_col))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filtered),
  colData = col_data,
  design = design_formula
)

# Set Reference Level
dds[[condition_col]] <- relevel(dds[[condition_col]],
                                ref = snakemake@params[["baselevel"]])

dds <- DESeq(dds)

saveRDS(dds, file = snakemake@output[["rds"]])