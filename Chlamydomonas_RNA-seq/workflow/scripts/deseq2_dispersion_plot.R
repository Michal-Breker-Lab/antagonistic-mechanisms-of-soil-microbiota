log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("ggplot2")

# Load data
dds <- readRDS(snakemake@input[["rds"]])

svg(snakemake@output[["dispersion_plot"]], width = 7, height = 7)
plotDispEsts(dds)
dev.off()

# Extract values
disp_values <- mcols(dds)$dispGeneEst
min_disp <- min(disp_values, na.rm = TRUE)
num_at_floor <- sum(disp_values == min_disp, na.rm = TRUE)
total_genes <- sum(!is.na(disp_values))

# Create a summary table
disp_summary <- data.frame(
  Metric = c("Min_Dispersion_Value",
             "Genes_at_Floor",
             "Total_Filtered_Genes",
             "Percent_at_Floor"),
  Value = c(min_disp,
            num_at_floor,
            total_genes,
            paste0(round((num_at_floor / total_genes) * 100, 2), "%"))
)

write.table(disp_summary,
            file = snakemake@output[["dispersion_stats"]],
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)