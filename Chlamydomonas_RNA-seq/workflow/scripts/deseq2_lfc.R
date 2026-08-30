log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(DESeq2)
library(ashr)

dds <- readRDS(snakemake@input[["rds"]])

condition_col <- snakemake@params[["condition"]]
baselevel     <- snakemake@params[["baselevel"]]

timepoints <- levels(dds[[condition_col]])
timepoints <- timepoints[timepoints != baselevel]

message("Running comparisons vs ", baselevel, " for: ", paste(timepoints, collapse = ", "))

results_list <- lapply(timepoints, function(tp) {
    message("Processing: ", tp, " vs ", baselevel)

    res <- results(dds, contrast = c(condition_col, tp, baselevel))

    res_shrunk <- lfcShrink(dds,
                            contrast = c(condition_col, tp, baselevel),
                            type     = "ashr",
                            res      = res)

    df <- as.data.frame(res_shrunk)
    df$gene       <- rownames(df)
    df$comparison <- paste0("T", tp, "_vs_T", baselevel)
    df
})

results_all <- do.call(rbind, results_list)

# save long format
write.table(results_all,
            snakemake@output[["long"]],
            sep       = "\t",
            row.names = FALSE,
            quote     = FALSE)

message("Done.")