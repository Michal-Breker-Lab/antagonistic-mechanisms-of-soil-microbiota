log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(ImpulseDE2)
library(dplyr)

# Use deseq2 RDS object to extract passed genes
dds <- readRDS(snakemake@input[["deseq2_rds"]])
genes_dds <- rownames(dds)

counts <- read.table(snakemake@input[["counts"]],
                     header = TRUE,
                     row.names = "gene",
                     check.names = FALSE)
counts <- as.matrix(round(counts))

anno <- read.table(snakemake@config[["samples"]],
                   sep = ",",
                   header = TRUE,
                   check.names = FALSE)

# Keep Genes that passed filtration during deseq2 step
counts_filtered <- counts[
  rownames(counts) %in% genes_dds,
]

# Log results
message(paste("Original genes:", nrow(counts)))
message(paste("Genes after strict timepoint filtering:", nrow(counts_filtered)))

time_col <- snakemake@params[["condition"]]
df_anno <- anno %>%
  transmute(
    Sample = sample_name,
    Condition = "case",
    Time = as.numeric(as.character(.data[[time_col]]))
  )

rownames(df_anno) <- df_anno$Sample

# Run ImpulseDE2
impulse_results <- runImpulseDE2(
  matCountData = counts_filtered,
  dfAnnotation = df_anno,
  boolCaseCtrl = FALSE,
  vecConfounders = NULL,
  scaNProc = snakemake@threads,
  boolIdentifyTransients = TRUE,
  scaQThres = as.numeric(snakemake@params[["scaQThres"]])
)

# Save Outputs
saveRDS(impulse_results,
        file = snakemake@output[["rds"]])
