#!/usr/bin/env Rscript
# Stage 11 - phylogenetically corrected tests of c3 presence and host association.
#
# WHY NOT A CHI-SQUARE: c3 presence is a binary trait distributed on a tree.
# Genomes are not independent observations. A single ancient loss inside a
# heavily sequenced clade would be counted as hundreds of independent events and
# return a spectacular, meaningless p-value. Everything here corrects for that.

suppressPackageStartupMessages({
  W <- Sys.getenv("W", "/mnt/LargeStorageNoBackup/Moshea/burkholderia_c3")
  .libPaths(c(file.path(W, "Rlib"), .libPaths()))
  library(ape); library(caper); library(phytools); library(phylolm)
})

RES <- file.path(W, "results")
TREE <- file.path(W, "trees", "chr1_core.treefile")
OUT <- file.path(RES, "phylo_stats.txt")
sink(OUT, split = TRUE)

cat("=== Stage 11: phylogenetic tests ===\n")
tr <- read.tree(TREE)
cat("tree tips:", Ntip(tr), "\n")

# ---- traits ----
cl <- read.delim(file.path(RES, "secondary_replicon_clusters.tsv"),
                 stringsAsFactors = FALSE)
c3 <- tapply(cl$is_c3 == "True", cl$accession, any)
hosts <- read.delim(file.path(RES, "host_categories.tsv"), stringsAsFactors = FALSE)

d <- data.frame(accession = tr$tip.label, stringsAsFactors = FALSE)
d$c3 <- as.integer(ifelse(is.na(c3[d$accession]), FALSE, c3[d$accession]))
d$host <- hosts$host_category[match(d$accession, hosts$accession)]
d$host[is.na(d$host)] <- "unknown"
d$organism <- hosts$organism_name[match(d$accession, hosts$accession)]

cat("\nc3 present:", sum(d$c3), " absent:", sum(!d$c3), "\n")
cat("\nhost categories on the tree:\n"); print(sort(table(d$host), decreasing = TRUE))

# tree must be dichotomous & rooted for several methods
# caper::comparative.data rejects a tree whose node labels repeat (IQ-TREE
# writes UFBoot values as node labels, so "100" appears hundreds of times).
# The support values are not used by any test here, so drop them.
tr$node.label <- NULL
tr <- multi2di(tr); tr$edge.length[tr$edge.length <= 0] <- 1e-8

# ---------------------------------------------------------------------------
# 1. Phylogenetic signal in c3 presence (Fritz & Purvis D)
#    D ~ 0  -> clumped as expected under Brownian evolution (clade-structured)
#    D ~ 1  -> randomly scattered across the tree (repeated independent loss)
#    D < 0  -> even more clumped than Brownian
# ---------------------------------------------------------------------------
cat("\n\n=== 1. PHYLOGENETIC SIGNAL IN c3 PRESENCE (Fritz & Purvis D) ===\n")
cd <- comparative.data(tr, d, names.col = "accession", na.omit = FALSE)
pd <- phylo.d(cd, binvar = c3, permut = 1000)
print(pd)
cat("\nInterpretation: D =", round(pd$DEstimate, 3),
    "\n  P(D = 0, Brownian/clumped) =", signif(pd$Pval0, 3),
    "\n  P(D = 1, random)           =", signif(pd$Pval1, 3), "\n")

# ---------------------------------------------------------------------------
# 2. How many independent gains/losses? Stochastic character mapping.
# ---------------------------------------------------------------------------
cat("\n\n=== 2. ANCESTRAL STATE RECONSTRUCTION (independent events) ===\n")
x <- setNames(factor(ifelse(d$c3 == 1, "present", "absent")), d$accession)
x <- x[tr$tip.label]
fit_er <- fitMk(tr, x, model = "ER")
fit_ard <- fitMk(tr, x, model = "ARD")
cat("ER  logLik:", round(logLik(fit_er), 2), " AIC:", round(AIC(fit_er), 2), "\n")
cat("ARD logLik:", round(logLik(fit_ard), 2), " AIC:", round(AIC(fit_ard), 2), "\n")
best <- if (AIC(fit_ard) < AIC(fit_er) - 2) "ARD" else "ER"
cat("model selected:", best, "\n")
sm <- make.simmap(tr, x, model = best, nsim = 100, message = FALSE)
ct <- describe.simmap(sm)
cat("\nmean transitions over 100 stochastic maps:\n"); print(round(ct$count[1, ], 2))
cat("\nThis is the number of INDEPENDENT gain/loss events, not the number of\n",
    "genomes lacking c3 -- the quantity that actually matters for asking whether\n",
    "loss is a recurrent, tolerated event.\n")

# ---------------------------------------------------------------------------
# 3. Host association, corrected for phylogeny
# ---------------------------------------------------------------------------
cat("\n\n=== 3. HOST ASSOCIATION (phylogenetic logistic regression) ===\n")
sub <- d[d$host != "unknown", ]
keep <- names(which(table(sub$host) >= 10))
sub <- sub[sub$host %in% keep, ]
cat("genomes used:", nrow(sub), " categories:", paste(keep, collapse = ", "), "\n")

if (nrow(sub) > 30 && length(keep) > 1) {
  tr2 <- drop.tip(tr, setdiff(tr$tip.label, sub$accession))
  sub <- sub[match(tr2$tip.label, sub$accession), ]
  rownames(sub) <- sub$accession
  ref <- if ("soil" %in% keep) "soil" else keep[1]
  sub$host <- relevel(factor(sub$host), ref = ref)
  cat("reference category:", ref, "\n")

  cat("\n--- naive (WRONG, shown only for contrast) ---\n")
  print(summary(glm(c3 ~ host, data = sub, family = binomial))$coefficients)
  cat("\n--- raw contingency table ---\n")
  print(table(sub$host, ifelse(sub$c3 == 1, "c3_present", "c3_absent")))

  cat("\n--- phylogenetic logistic regression (Ives & Garland) ---\n")
  fit <- tryCatch(
    phyloglm(c3 ~ host, phy = tr2, data = sub, method = "logistic_MPLE",
             btol = 30, log.alpha.bound = 6),
    error = function(e) { cat("phyloglm failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(fit)) {
    print(summary(fit))
    co <- summary(fit)$coefficients
    cat("\nodds ratios with 95% CI:\n")
    for (i in seq_len(nrow(co))) {
      est <- co[i, 1]; se <- co[i, 2]
      cat(sprintf("  %-28s OR=%6.3f  [%6.3f, %6.3f]  p=%.4g\n",
                  rownames(co)[i], exp(est), exp(est - 1.96 * se),
                  exp(est + 1.96 * se), co[i, 4]))
    }
  }
} else {
  cat("insufficient data for host regression\n")
}

# ---------------------------------------------------------------------------
# 3b. Species tree vs c3 tree: vertical inheritance or horizontal movement?
# ---------------------------------------------------------------------------
cat("\n\n=== 3b. TREE CONGRUENCE (chr1 species tree vs c3 tree) ===\n")
C3TREE <- file.path(W, "trees", "c3_core.treefile")
if (file.exists(C3TREE)) {
  t2 <- read.tree(C3TREE)
  shared <- intersect(tr$tip.label, t2$tip.label)
  cat("shared tips:", length(shared), "\n")
  if (length(shared) >= 5) {
    a <- drop.tip(tr, setdiff(tr$tip.label, shared))
    b <- drop.tip(t2, setdiff(t2$tip.label, shared))
    rf <- dist.topo(a, b, method = "PH85")
    maxrf <- 2 * (length(shared) - 3)
    cat("Robinson-Foulds distance:", rf, "of a maximum", maxrf,
        sprintf("(%.1f%%)\n", 100 * rf / maxrf))
    cat("normalised congruence:", sprintf("%.3f\n", 1 - rf / maxrf))

    # null: RF between the species tree and random trees on the same tips
    null <- replicate(200, dist.topo(a, rtree(length(shared), tip.label = shared),
                                     method = "PH85"))
    cat("RF vs random trees: mean", sprintf("%.1f", mean(null)),
        " sd", sprintf("%.1f", sd(null)), "\n")
    z <- (rf - mean(null)) / sd(null)
    cat("z =", sprintf("%.2f", z), "\n")
    cat(if (rf < mean(null) - 3 * sd(null))
      "  -> trees are significantly MORE congruent than random: c3 largely\n     tracks the chromosome, i.e. vertical inheritance dominates.\n"
      else
      "  -> congruence is not distinguishable from random: c3 history is\n     decoupled from the chromosome, i.e. horizontal movement dominates.\n")
  }
} else {
  cat("c3 tree not available - skipping congruence test\n")
}

# ---------------------------------------------------------------------------
# 4. The sampling confound, stated numerically
# ---------------------------------------------------------------------------
cat("\n\n=== 4. SAMPLING CONFOUND ===\n")
tb <- sort(table(d$organism), decreasing = TRUE)
cat("top 10 organism labels on the tree:\n"); print(head(tb, 10))
cat("\ntop-10 share of tips:", round(100 * sum(head(tb, 10)) / Ntip(tr), 1), "%\n")
cat("Sequencing effort in this genus is driven by clinical and biodefense\n",
    "priorities, not by random sampling of Burkholderia diversity. Any host\n",
    "association must be read with that in mind; dereplication and the\n",
    "phylogenetic correction above mitigate it but do not remove it.\n")

sink()
cat("wrote", OUT, "\n")
cat("STAGE11_DONE\n")
