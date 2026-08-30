log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")
library("ggplot2")

# Load data
dds <- readRDS(snakemake@input[["rds"]])

# Selection of Transformation Method
method <- snakemake@params[["transformation"]]

if (method == "rlog") {
  transformed_data <- rlog(dds, blind = TRUE)
} else if (method == "vst") {
  transformed_data <- vst(dds, blind = TRUE)
} else {
  stop(paste0("Invalid transformation method: '", method,
              "'. Please use 'vst' or 'rlog' in your config file."))
}

# PCA Calculation
ntop_val <- as.numeric(snakemake@params[["ntop"]])

pcaData <- plotPCA(
  transformed_data,
  intgroup = snakemake@params[["condition"]],
  ntop = ntop_val,
  returnData = TRUE
)

# 3. Plotting
percentVar <- round(100 * attr(pcaData, "percentVar"))
condition_col <- snakemake@params[["condition"]]

p <- ggplot(pcaData, aes(PC1, PC2, color = .data[[condition_col]])) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  labs(
    color = condition_col,
    #title = "PCA of normalized counts",
  )

# Save output
ggsave(
  snakemake@output[["pca_plot"]],
  p,
  width = 6,
  height = 5,
  device = "svg",
  bg = "white"
)