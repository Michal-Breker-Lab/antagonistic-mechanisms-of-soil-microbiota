log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(Mfuzz)

exps_up <- read.table(snakemake@input[["zscores"]],
                      header = TRUE, row.names = 1,
                      sep = "\t", check.names = FALSE)

expr_mat <- as.matrix(exps_up)
storage.mode(expr_mat) <- "numeric"

eset <- ExpressionSet(assayData = expr_mat)

set.seed(as.numeric(snakemake@params[["seed"]]))
m <- mestimate(eset)
nclusters <- as.numeric(snakemake@params[["Nclusters"]])
cl <- mfuzz(eset, c = nclusters, m = m)

centers <- cl$centers
write.table(centers, file = snakemake@output[["centers"]],
            sep = "\t", quote = FALSE, col.names = NA)

membership <- cl$membership
write.table(membership, file = snakemake@output[["membership"]],
            sep = "\t", quote = FALSE, col.names = NA)

clusters <- data.frame(Gene = names(cl$cluster), Cluster = cl$cluster)
write.table(clusters, file = snakemake@output[["clusters"]],
            sep = "\t", quote = FALSE, row.names = FALSE)