log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(ImpulseDE2)
library(dplyr)

impulse_rds <- readRDS(snakemake@input[["impulsede2_rds"]])

# Per-gene fit results: p, padj (the q-value), isSignificant (padj < scaQThres),
# transient/monotonous shape flags, ...
results <- as.data.frame(impulse_rds$dfImpulseDE2Results)

# Re-derive the transient/transition up/down classification exactly as the
# downstream zscores subsetting does (impulseDE2_zscores.R), so this table
# reflects the gene groups that are actually used later in the pipeline.
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = impulse_rds,
  strCondition           = "case",
  boolIdentifyTransients = TRUE,
  scaQThres              = as.numeric(snakemake@params[["scaQThres"]])
)

df_groups <- stack(lsHeatmaps$lsvecGeneGroups)
colnames(df_groups) <- c("Gene", "category")
# category is e.g. transient_up / transition_up / transient_down / transition_down
df_groups$direction <- sub(".*_", "", df_groups$category)

results <- results %>%
  left_join(df_groups, by = "Gene")

# Append annotation columns, same set as add_annotations_to_table.py.
annotation_cols <- c("geneSymbol", "Description", "Predalgo", "GO",
                     "Best-hit-arabi-name")
annotations <- read.delim(snakemake@input[["all_annotations"]],
                          check.names = FALSE, quote = "",
                          comment.char = "", row.names = 1)
annotations$Gene <- rownames(annotations)
annotations <- annotations[, c("Gene", annotation_cols)]

results <- results %>%
  left_join(annotations, by = "Gene")

write.table(results,
            file = snakemake@output[["table"]],
            sep = "\t", quote = FALSE, row.names = FALSE)
