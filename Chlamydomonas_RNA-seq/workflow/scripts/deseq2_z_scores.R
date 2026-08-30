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
  transformed_data <- rlog(dds, blind = FALSE)
} else if (method == "vst") {
  transformed_data <- vst(dds, blind = FALSE)
} else {
  stop(paste0("Invalid transformation method: '", method,
              "'. Please use 'vst' or 'rlog' in your config file."))
}

mat <- assay(transformed_data)

metadata <- as.data.frame(colData(dds))
condition_col <- snakemake@params[["condition"]]

# Collapse samples into condition means
mean_mat <- t(apply(mat, 1, function(x) {
  tapply(x, metadata[[condition_col]], mean)
}))

# Calculate Z-score per condition
exps <- t(scale(t(mean_mat)))

# Remove NaNs
exps <- exps[complete.cases(exps), ]

exps_df <- as.data.frame(exps)
# Save as a TSV file
write.table(exps_df,
            file = snakemake@output[["zscores"]],
            sep = "\t", 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = NA) # This keeps the 'Timepoint' headers aligned