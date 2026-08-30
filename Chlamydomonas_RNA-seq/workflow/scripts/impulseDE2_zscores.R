log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(ImpulseDE2)
library(dplyr)

direction <- snakemake@params[["direction"]]

zscores <- read.table(snakemake@input[["zscores"]],
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE)

impulse_rds <- readRDS(snakemake@input[["impulsede2_rds"]])

lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = impulse_rds,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = as.numeric(snakemake@params[["scaQThres"]])
)

df_groups <- stack(lsHeatmaps$lsvecGeneGroups)

colnames(df_groups) <- c("gene", "category")

# 2. Define the genes for the 'up' group
# We include both 'transition_up' and 'transient_up'
sub_gene_ids <- df_groups %>%
  filter(category %in% c(paste0("transient_", direction),
                         paste0("transition_", direction))) %>%
  pull(gene)

zscores_sub <- zscores[rownames(zscores) %in% sub_gene_ids, ]

write.table(zscores_sub,
            file = snakemake@output[["zscores_direction"]],
            sep = "\t", quote = FALSE, col.names = NA)
