#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(limma)})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output"); sink(log, type = "message")

in_tsv <- snakemake@input[["matrix"]]
outdir <- dirname(snakemake@output[[1]])
gff    <- snakemake@input[["gff"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

MIN_VALID <- as.integer(snakemake@params[["min_valid"]])
MIN_ON    <- as.integer(snakemake@params[["min_on"]])
FDR       <- as.numeric(snakemake@params[["fdr"]])
LFC       <- as.numeric(snakemake@params[["lfc"]])
DROP      <- as.character(snakemake@params[["drop"]])
DAY       <- as.integer(snakemake@params[["day"]])
N_RUNS    <- as.integer(snakemake@params[["n_runs"]])
N_MF6     <- as.integer(snakemake@params[["n_mf6"]])
MIN_CORE  <- as.integer(snakemake@params[["min_core"]])

ABSENT_IN_MUTANTS <- snakemake@params[["absent"]]
contig_of <- local({
  g <- read.delim(gff, header = FALSE,
                  comment.char = "#", stringsAsFactors = FALSE, quote = "")
  g <- g[g[[3]] == "CDS", ]
  setNames(g[[1]], sub(".*locus_tag=([^;]+).*", "\\1", g[[9]]))
})

x <- read.delim(in_tsv, check.names = FALSE)
stopifnot("Protein" %in% names(x))
ann_cols <- intersect(c("Protein", "Gene", "Protein Length", "Combined Total Peptides",
                        "Combined Unique Spectral Count", "Protein Probability",
                        "Description"), names(x))

all_d2 <- setdiff(grep(sprintf("_d%d_[1-4]$", DAY), names(x), value = TRUE), DROP)
mf6_d2 <- grep(sprintf("^MF6_(alone|withC)_d%d_[1-4]$", DAY), names(x), value = TRUE)
stopifnot(length(all_d2) == N_RUNS, length(mf6_d2) == N_MF6)

A <- as.matrix(x[, all_d2]); storage.mode(A) <- "double"
A[A <= 0] <- NA
A <- log2(A)
rownames(A) <- x$Protein

cat(sprintf("input      : %s\n", in_tsv))
cat(sprintf("excluded   : %s (failed injection, QC_report.md)\n", DROP))
cat(sprintf("proteins   : %d\nday-2 runs : %d (core/offsets)   MF6 runs: %d (fit)\n",
            nrow(A), ncol(A), length(mf6_d2)))

core <- rowSums(is.na(A)) == 0 &
        !(contig_of[rownames(A)] %in% ABSENT_IN_MUTANTS)
cat(sprintf("CORE       : %d of %d proteins with a real value in all %d day-2 runs (%s excluded)\n",
            sum(core), nrow(A), ncol(A), paste(ABSENT_IN_MUTANTS, collapse = "/")))
stopifnot(sum(core) > MIN_CORE)

offset <- apply(A[core, , drop = FALSE], 2, median)
A <- sweep(A, 2, offset, "-") + median(offset)
cat(sprintf("centring   : core-set medians %.2f-%.2f (spread %.2f log2) -> %.2f (all)\n",
            min(offset), max(offset), max(offset) - min(offset), median(offset)))

E <- A[, mf6_d2, drop = FALSE]

group <- factor(sub("_[1-4]$", "", colnames(E)),
                levels = c("MF6_alone_d2", "MF6_withC_d2"))
nvalid <- sapply(levels(group), function(g)
  rowSums(!is.na(E[, group == g, drop = FALSE])))
colnames(nvalid) <- levels(group)
quantified <- nvalid >= MIN_VALID
keep <- rowSums(quantified) >= 1
cat(sprintf("kept       : %d proteins quantified in >=%d/4 replicates of >=1 group\n",
            sum(keep), MIN_VALID))

Ek <- E[keep, , drop = FALSE]; qk <- quantified[keep, , drop = FALSE]
nvk <- nvalid[keep, , drop = FALSE]; annk <- x[keep, ann_cols, drop = FALSE]

design <- model.matrix(~ 0 + group); colnames(design) <- levels(group)
cm <- makeContrasts(coculture_d2 = MF6_withC_d2 - MF6_alone_d2, levels = design)
fitable <- qk[, "MF6_alone_d2"] & qk[, "MF6_withC_d2"]
fit2 <- eBayes(contrasts.fit(lmFit(Ek[fitable, , drop = FALSE], design), cm),
               trend = TRUE, robust = TRUE)
cat(sprintf("limma      : %d residual df, prior df %.1f\n",
            fit2$df.residual[1], fit2$df.prior[1]))

grp_mean <- sapply(levels(group), function(g)
  rowMeans(Ek[, group == g, drop = FALSE], na.rm = TRUE))
colnames(grp_mean) <- levels(group)

tt <- topTable(fit2, coef = "coculture_d2", number = Inf, sort.by = "none")
res <- data.frame(annk[fitable, , drop = FALSE],
                  logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                  P.Value = tt$P.Value,
                  adj.P.Val = p.adjust(tt$P.Value, method = "BH"), B = tt$B,
                  group_A = "MF6_withC_d2", group_B = "MF6_alone_d2",
                  mean_A = grp_mean[fitable, "MF6_withC_d2"],
                  mean_B = grp_mean[fitable, "MF6_alone_d2"],
                  n_valid_A = nvk[fitable, "MF6_withC_d2"],
                  n_valid_B = nvk[fitable, "MF6_alone_d2"],
                  check.names = FALSE)
res$significant <- res$adj.P.Val < FDR
res$sig_and_lfc <- res$significant & abs(res$logFC) >= LFC
res$direction   <- ifelse(res$logFC > 0, "up", "down")
res <- res[order(res$adj.P.Val, res$P.Value), ]

write.table(res, file.path(outdir, "DE_coculture_d2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("\n  %-22s tested %5d   FDR<%.2f %5d   |logFC|>=%g %5d  (up %d / down %d)\n",
            "coculture_d2", nrow(res), FDR, sum(res$significant), LFC,
            sum(res$sig_and_lfc), sum(res$sig_and_lfc & res$logFC > 0),
            sum(res$sig_and_lfc & res$logFC < 0)))

summ <- data.frame(contrast = "coculture_d2", group_A = "MF6_withC_d2",
                   group_B = "MF6_alone_d2", tested = nrow(res),
                   sig_FDR05 = sum(res$significant),
                   sig_and_lfc1 = sum(res$sig_and_lfc),
                   up = sum(res$sig_and_lfc & res$logFC > 0),
                   down = sum(res$sig_and_lfc & res$logFC < 0))
write.table(summ, file.path(outdir, "DE_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pctile <- function(v, ref) {
  ref <- ref[is.finite(ref)]
  if (!length(ref)) return(NA_real_)
  100 * sapply(v, function(z) mean(ref < z, na.rm = TRUE))
}

onoff <- list()
for (side in 1:2) {
  on_g  <- levels(group)[c(2, 1)[side]]
  off_g <- levels(group)[c(1, 2)[side]]
  sel <- nvk[, on_g] >= MIN_ON & nvk[, off_g] == 0
  if (!any(sel)) next
  mu <- grp_mean[sel, on_g]
  onoff[[length(onoff) + 1]] <- data.frame(
    annk[sel, , drop = FALSE],
    present_in = sub("^MF6_", "", on_g), absent_from = sub("^MF6_", "", off_g),
    n_valid_present = nvk[sel, on_g], n_valid_absent = nvk[sel, off_g],
    mean_log2_LFQ_present = round(mu, 3),
    pct_in_own_group      = round(pctile(mu, grp_mean[, on_g]), 1),
    pct_in_opposite_group = round(pctile(mu, grp_mean[, off_g]), 1),
    check.names = FALSE)
}
onoff <- do.call(rbind, onoff)
if (!is.null(onoff)) {
  onoff$evidence <- cut(onoff$pct_in_opposite_group, c(-Inf, 20, 70, Inf),
                        labels = c("weak (near detection limit)", "moderate",
                                   "strong (well within range)"))
  onoff <- onoff[order(-onoff$pct_in_opposite_group), ]
  f <- file.path(outdir, "on_off_withC_d2_vs_alone_d2.tsv")
  write.table(onoff, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  ev <- table(factor(onoff$evidence, levels = levels(onoff$evidence)))
  cat(sprintf("  %-22s %5d rows  (strong %d / moderate %d / weak %d)\n",
              "on/off", nrow(onoff), ev[3], ev[2], ev[1]))
}

write.table(data.frame(annk, round(Ek, 4), check.names = FALSE),
            file.path(outdir, "MF6_log2_maxlfq_core_centred.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("\nwritten    -> %s\n", outdir))
