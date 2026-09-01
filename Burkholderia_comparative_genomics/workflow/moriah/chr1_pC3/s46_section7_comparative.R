#!/usr/bin/env Rscript
# Section 7 comparative analyses on the rebuilt 763-tip chromosome-1 tree.
#
# Reproduces the methods the report already declares, on the new tree:
#   1. Stochastic character mapping of pC3 presence (ARD vs ER chosen by AIC,
#      100 maps) -> counts of loss and regain events.
#   2. Fritz & Purvis D for pC3 presence, with both null distributions.
#   3. Phylogenetic logistic regression of pC3 presence on host category,
#      soil as the reference level, beside the uncorrected GLM so the two can
#      be compared -- the point of section 7.2 is that they disagree.
#
# The tree is UNROOTED (IQ-TREE writes a trifurcating root). simmap and D both
# need a rooted tree, so it is midpoint-rooted here; that is a computational
# requirement, and section 7 already states the rooting is a display convention
# rather than a claim. The clan and patristic tests elsewhere avoid it entirely.
suppressPackageStartupMessages({
  library(ape); library(phytools); library(caper); library(phylolm)
})
set.seed(20260825)

args <- commandArgs(trailingOnly = TRUE)
tree_f <- args[1]; calls_f <- args[2]; host_f <- args[3]; out_f <- args[4]

# Results are appended to the output file as each stage completes. The first
# attempt (job 45940924) ran make.simmap for ~70 minutes, then died in a later
# stage and wrote nothing at all.
LINES <- c("analysis\tterm\tvalue\tnote")
emit <- function(...) {
  LINES <<- c(LINES, sprintf(...))
  writeLines(LINES, out_f)
}

tr <- read.tree(tree_f)
cat(sprintf("tree: %d tips\n", Ntip(tr)))
tr <- midpoint.root(multi2di(tr))
tr$edge.length[tr$edge.length <= 0] <- 1e-8

calls <- read.delim(calls_f, stringsAsFactors = FALSE)
host  <- read.delim(host_f,  stringsAsFactors = FALSE)

pres <- setNames(as.integer(calls$c3_present == "True"), calls$accession)
pres <- pres[tr$tip.label]
stopifnot(!any(is.na(pres)))
cat(sprintf("pC3 present on tree: %d / %d\n", sum(pres), length(pres)))

# ---- 1. model choice, then stochastic character mapping --------------------
x <- setNames(factor(ifelse(pres == 1, "present", "absent")), tr$tip.label)
fit_er  <- fitMk(tr, x, model = "ER")
fit_ard <- fitMk(tr, x, model = "ARD")
aic <- c(ER = AIC(fit_er), ARD = AIC(fit_ard))
best <- names(which.min(aic))
cat(sprintf("AIC  ER %.1f   ARD %.1f   -> %s\n", aic["ER"], aic["ARD"], best))

emit("simmap\tmodel_chosen\t%s\tAIC ER %.1f vs ARD %.1f", best, aic["ER"], aic["ARD"])
cnt_f <- sub("\\.tsv$", "_simmap_counts.rds", out_f)
if (file.exists(cnt_f)) {
  cnt <- readRDS(cnt_f)          # ~70 min to recompute; reuse it
  cat("reusing saved simmap counts\n")
} else {
  maps <- make.simmap(tr, x, model = best, nsim = 100, message = FALSE)
  cnt <- t(sapply(maps, function(m) {
    ct <- countSimmap(m)$Tr
    c(loss = ct["present", "absent"], gain = ct["absent", "present"])
  }))
}
mean_loss <- mean(cnt[, "loss"]); mean_gain <- mean(cnt[, "gain"])
cat(sprintf("simmap over 100 maps: present->absent %.1f (sd %.1f), absent->present %.1f (sd %.1f)\n",
            mean_loss, sd(cnt[, "loss"]), mean_gain, sd(cnt[, "gain"])))
emit("simmap\tlosses_present_to_absent\t%.1f\tmean over 100 maps, sd %.1f", mean_loss, sd(cnt[,"loss"]))
emit("simmap\tgains_absent_to_present\t%.1f\tmean over 100 maps, sd %.1f", mean_gain, sd(cnt[,"gain"]))
saveRDS(cnt, sub("\\.tsv$", "_simmap_counts.rds", out_f))

# ---- 2. Fritz & Purvis D ---------------------------------------------------
# IQ-TREE stores UFBoot support as internal node labels, and many nodes share
# the value 100; caper reads those as duplicated labels and refuses the tree
# ("Labels duplicated between tips and nodes in phylogeny"). The support values
# play no part in these analyses, so drop them.
tr_nl <- tr
tr_nl$node.label <- NULL
dd <- data.frame(tip = tr_nl$tip.label, pc3 = as.integer(pres))
cd <- comparative.data(tr_nl, dd, names.col = "tip")
d <- tryCatch(phylo.d(cd, binvar = pc3, permut = 1000), error = function(e) {
  emit("fritz_purvis\tD\tNA\tFAILED: %s", conditionMessage(e)); NULL })
if (!is.null(d)) {
  cat(sprintf("Fritz & Purvis D = %.3f   P(random) = %.3f   P(Brownian) = %.3f\n",
              d$DEstimate, d$Pval1, d$Pval0))
  emit("fritz_purvis\tD\t%.3f\tP(random)=%.3f P(Brownian)=%.3f", d$DEstimate, d$Pval1, d$Pval0)
}

# ---- 3. host association, corrected and uncorrected ------------------------
h <- setNames(host$host_category, host$accession)[tr$tip.label]
keep <- !is.na(h) & h %in% names(which(table(h) >= 20)) & h != "unknown"
cat(sprintf("host-modelled tips: %d in %d categories\n", sum(keep), length(unique(h[keep]))))
dh <- data.frame(pc3 = as.integer(pres[keep]),
                 host = relevel(factor(h[keep]), ref = "soil"))
rownames(dh) <- tr$tip.label[keep]
trh <- keep.tip(tr_nl, tr_nl$tip.label[keep])

emit("host\tn_tips_modelled\t%d\tcategories with n>=20, unknown excluded", sum(keep))
glm_fit <- glm(pc3 ~ host, data = dh, family = binomial)
sg <- summary(glm_fit)$coefficients
for (nm in rownames(sg))
  emit("glm_uncorrected\t%s\t%.3f\tp=%.3g", nm, sg[nm, 1], sg[nm, 4])

pg <- tryCatch(phyloglm(pc3 ~ host, data = dh, phy = trh,
                        method = "logistic_MPLE", btol = 50),
               error = function(e) { emit("phyloglm\t(all)\tNA\tFAILED: %s",
                                          conditionMessage(e)); NULL })
if (is.null(pg)) { cat("phyloglm failed; partial results written\n"); quit(status = 0) }
sp <- summary(pg)$coefficients

for (nm in rownames(sp))
  emit("phyloglm\t%s\t%.3f\tp=%.3g", nm, sp[nm, 1], sp[nm, ncol(sp)])
cat(sprintf("\nwrote %s\n", out_f))
