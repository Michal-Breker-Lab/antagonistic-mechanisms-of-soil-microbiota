#!/usr/bin/env Rscript

suppressPackageStartupMessages({library(limma)})

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output"); sink(log, type = "message")

in_tsv  <- snakemake@input[["matrix"]]
de_root <- dirname(snakemake@output[[1]])
gff     <- snakemake@input[["gff"]]

DROP <- as.character(snakemake@params[["drop"]])
FDR  <- as.numeric(snakemake@params[["fdr"]])
LFC  <- as.numeric(snakemake@params[["lfc"]])
DAY  <- as.integer(snakemake@params[["day"]])
N_RUNS   <- as.integer(snakemake@params[["n_runs"]])
MIN_CORE <- as.integer(snakemake@params[["min_core"]])

ABSENT_IN_MUTANTS <- snakemake@params[["absent"]]
contig_of <- local({
  g <- read.delim(gff, header = FALSE, comment.char = "#",
                  stringsAsFactors = FALSE, quote = "")
  g <- g[g[[3]] == "CDS", ]
  lt <- sub(".*locus_tag=([^;]+).*", "\\1", g[[9]])
  setNames(g[[1]], lt)
})

x <- read.delim(in_tsv, check.names = FALSE, quote = "")
ann_cols <- intersect(c("Protein", "Gene", "Protein Length",
                        "Combined Total Peptides", "Combined Unique Spectral Count",
                        "Protein Probability", "Description"), colnames(x))

samples <- grep("^(MF6|27D6|34F7)_(alone|withC)_d2_[1-4]$", colnames(x), value = TRUE)
samples <- setdiff(samples, DROP)
stopifnot(length(samples) == N_RUNS)

E <- as.matrix(x[, samples])
storage.mode(E) <- "double"
E[E <= 0] <- NA
E <- log2(E)
rownames(E) <- x$Protein

sf <- function(g) sub("^(?=[0-9])", "S", g, perl = TRUE)
un <- function(g) sub("^S(?=[0-9])", "", g, perl = TRUE)

group <- factor(sf(sub("_[1-4]$", "", samples)),
                levels = sf(c("MF6_alone_d2", "27D6_alone_d2", "34F7_alone_d2",
                              "MF6_withC_d2", "27D6_withC_d2", "34F7_withC_d2")))
nrep  <- table(group)

cat(sprintf("input      : %s\n", in_tsv))
cat(sprintf("excluded   : %s (failed injection, QC_report.md)\n", DROP))
cat(sprintf("proteins   : %d\nsamples    : %d (%s)\n", nrow(E), ncol(E),
            paste(sprintf("%s=%d", un(levels(group)), nrep), collapse = ", ")))

core <- rowSums(is.na(E)) == 0 &
        !(contig_of[rownames(E)] %in% ABSENT_IN_MUTANTS)
cat(sprintf("CORE       : %d of %d proteins with a real value in all %d samples (%s excluded)\n",
            sum(core), nrow(E), ncol(E), paste(ABSENT_IN_MUTANTS, collapse = "/")))
stopifnot(sum(core) > MIN_CORE)

offset <- apply(E[core, , drop = FALSE], 2, median)
E <- sweep(E, 2, offset, "-") + median(offset)
cat(sprintf("centring   : core-set medians %.2f-%.2f (spread %.2f log2) -> %.2f (all)\n",
            min(offset), max(offset), max(offset) - min(offset), median(offset)))

core_tsv <- file.path(de_root, "day2_core_proteins.tsv")
dir.create(de_root, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(x[core, ann_cols, drop = FALSE],
                       mean_log2_LFQ_centred = round(rowMeans(E[core, , drop = FALSE]), 3),
                       check.names = FALSE),
            core_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

write.table(data.frame(sample = samples, group = un(as.character(group)),
                       median_offset = round(offset, 4),
                       n_quantified = colSums(!is.na(E))),
            file.path(de_root, "day2_median_offsets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

nvalid <- sapply(levels(group), function(g)
  rowSums(!is.na(E[, group == g, drop = FALSE])))
colnames(nvalid) <- levels(group)
min_valid <- rep(2L, length(nrep))
MIN_ON    <- as.integer(snakemake@params[["min_on"]])
names(min_valid) <- levels(group)
quantified <- sweep(nvalid, 2, min_valid, ">=")

cat("min valid  : ",
    paste(sprintf("%s>=%d", un(levels(group)), min_valid), collapse = ", "), "\n", sep = "")
cat(sprintf("quantified : %s\n",
            paste(sprintf("%s=%d", un(levels(group)), colSums(quantified)), collapse = ", ")))

keep <- rowSums(quantified) >= 1
cat(sprintf("kept       : %d proteins quantified in >=1 group\n", sum(keep)))

Ek   <- E[keep, , drop = FALSE]
qk   <- quantified[keep, , drop = FALSE]
nvk  <- nvalid[keep, , drop = FALSE]
annk <- x[keep, ann_cols, drop = FALSE]

absent <- contig_of[annk$Protein] %in% ABSENT_IN_MUTANTS
absent[is.na(absent)] <- FALSE
cat(sprintf("censored   : %d proteins on %s - absent from both mutant genomes, "
            , sum(absent), paste(ABSENT_IN_MUTANTS, collapse = "/")))
cat("excluded from every contrast below\n")

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit  <- lmFit(Ek, design)

grp_mean <- sapply(levels(group), function(g)
  rowMeans(Ek[, group == g, drop = FALSE], na.rm = TRUE))
colnames(grp_mean) <- levels(group)

STRAINS <- as.character(snakemake@params[["mutants"]])

CONTRASTS <- function(m) list(
  c(sprintf("DE_%s_withC_vs_alone", m),     sf(sprintf("%s_withC_d2", m)), sf(sprintf("%s_alone_d2", m))),
  c(sprintf("DE_%s_withC_vs_MF6_withC", m), sf(sprintf("%s_withC_d2", m)), "MF6_withC_d2"),
  c(sprintf("DE_%s_alone_vs_MF6_alone", m), sf(sprintf("%s_alone_d2", m)), "MF6_alone_d2"))

pairs  <- do.call(c, lapply(STRAINS, CONTRASTS))
names(pairs) <- vapply(pairs, `[`, "", 1)

cm <- makeContrasts(contrasts = vapply(pairs, function(z) paste(z[2], "-", z[3]), ""),
                    levels = design)
colnames(cm) <- names(pairs)
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)
cat(sprintf("limma      : %d residual df, prior df %.1f\n",
            fit$df.residual[1], fit2$df.prior[1]))

summary_rows <- list()
onoff_rows   <- list()
for (m in STRAINS) {
  outdir <- file.path(de_root, m)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n%s -> %s\n", m, outdir))

  for (z in CONTRASTS(m)) {
    cn <- z[1]; A <- z[2]; B <- z[3]
    tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    ok <- qk[, A] & qk[, B] & !is.na(tt$logFC) & !absent

    res <- data.frame(annk,
                      logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                      P.Value = tt$P.Value, adj.P.Val = NA_real_, B = tt$B,
                      group_A = un(A), group_B = un(B),
                      mean_A = grp_mean[, A], mean_B = grp_mean[, B],
                      n_valid_A = nvk[, A], n_valid_B = nvk[, B],
                      check.names = FALSE)
    res$adj.P.Val[ok] <- p.adjust(res$P.Value[ok], method = "BH")
    res <- res[ok, , drop = FALSE]
    res$significant <- res$adj.P.Val < FDR
    res$sig_and_lfc <- res$significant & abs(res$logFC) >= LFC
    res$direction   <- ifelse(res$logFC > 0, "up", "down")
    res <- res[order(res$adj.P.Val, res$P.Value), ]

    write.table(res, file.path(outdir, sprintf("%s.tsv", cn)),
                sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    summary_rows[[cn]] <- data.frame(
      strain = m, contrast = cn, group_A = un(A), group_B = un(B),
      tested = nrow(res),
      sig_FDR05 = sum(res$significant, na.rm = TRUE),
      sig_and_lfc1 = sum(res$sig_and_lfc, na.rm = TRUE),
      up = sum(res$sig_and_lfc & res$logFC > 0, na.rm = TRUE),
      down = sum(res$sig_and_lfc & res$logFC < 0, na.rm = TRUE))
    s <- summary_rows[[cn]]
    cat(sprintf("  %-30s tested %5d   FDR<%.2f %5d   |logFC|>=%g %5d  (up %d / down %d)\n",
                sub("^DE_", "", cn), s$tested, FDR, s$sig_FDR05, LFC,
                s$sig_and_lfc1, s$up, s$down))
  }

  for (z in CONTRASTS(m)) {
    cn <- z[1]; A <- z[2]; B <- z[3]
    for (side in 1:2) {
      on_g  <- if (side == 1) A else B
      off_g <- if (side == 1) B else A
      sel   <- nvk[, on_g] >= MIN_ON & nvk[, off_g] == 0 & !absent
      if (!any(sel)) next
      mu  <- grp_mean[sel, on_g]
      ref <- grp_mean[, off_g]; ref <- ref[is.finite(ref)]
      pct <- 100 * sapply(mu, function(z_) mean(ref < z_, na.rm = TRUE))
      onoff_rows[[length(onoff_rows) + 1]] <- data.frame(
        annk[sel, , drop = FALSE], contrast = cn, outdir = m,
        present_in = un(on_g), absent_from = un(off_g),
        n_valid_present = nvk[sel, on_g], n_valid_absent = nvk[sel, off_g],
        mean_log2_LFQ_present = round(mu, 3),
        pct_in_own_group      = round(100 * sapply(mu, function(z_)
                                  mean(grp_mean[, on_g] < z_, na.rm = TRUE)), 1),
        pct_in_opposite_group = round(pct, 1),
        check.names = FALSE)
    }
  }

  sm <- do.call(rbind, summary_rows[vapply(CONTRASTS(m), `[`, "", 1)])
  write.table(sm, file.path(outdir, "DE_summary.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cols <- samples[group %in% c(sf(sprintf("%s_alone_d2", m)), sf(sprintf("%s_withC_d2", m)),
                               "MF6_alone_d2", "MF6_withC_d2")]
  write.table(data.frame(annk, round(Ek[, cols, drop = FALSE], 4), check.names = FALSE),
              file.path(outdir, sprintf("%s_log2_maxlfq_core_centred.tsv", m)),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

oo <- do.call(rbind, onoff_rows)
if (!is.null(oo)) {
  oo$evidence <- cut(oo$pct_in_opposite_group, c(-Inf, 20, 70, Inf),
                     labels = c("weak (near detection limit)", "moderate",
                                "strong (well within range)"))
  cat("\non/off:\n")
  for (cn in unique(oo$contrast)) {
    sub <- oo[oo$contrast == cn, , drop = FALSE]
    sub <- sub[order(-sub$pct_in_opposite_group), ]
    dest <- file.path(de_root, sub$outdir[1],
                      sprintf("on_off_%s.tsv", sub("^DE_", "", cn)))
    sub$contrast <- NULL; sub$outdir <- NULL
    write.table(sub, dest, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    ev <- table(factor(sub$evidence, levels = levels(oo$evidence)))
    cat(sprintf("  %-30s %4d rows  (strong %d / moderate %d / weak %d)  -> %s\n",
                sub("^DE_", "", cn), nrow(sub), ev[3], ev[2], ev[1], basename(dest)))
  }
}

write.table(data.frame(annk, round(Ek, 4), check.names = FALSE),
            file.path(de_root, "day2_log2_maxlfq_core_centred.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

write.table(do.call(rbind, summary_rows), file.path(de_root, "DE_summary_mutants_day2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\ncore set   -> %s\nmatrix     -> %s\nsummary    -> %s\n", core_tsv,
            file.path(de_root, "day2_log2_maxlfq_core_centred.tsv"),
            file.path(de_root, "DE_summary_mutants_day2.tsv")))
