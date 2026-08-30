#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output"); sink(log, type = "message")

DE_DIR    <- snakemake@params[["de_dir"]]
OUT_DIR   <- dirname(snakemake@output[[1]])
EGGNOG    <- snakemake@input[["eggnog"]]
KEGGNAME  <- snakemake@input[["kegg"]]
KEGGCLASS <- snakemake@input[["kegg_classes"]]
DROP_NONBACT  <- isTRUE(snakemake@params[["drop_nb"]])
DROP_OVERVIEW <- isTRUE(snakemake@params[["drop_ov"]])
DROP_RETIRED  <- isTRUE(snakemake@params[["drop_ret"]])

MIN_TERM_SIZE <- as.integer(snakemake@params[["min_term"]])
MIN_HITS      <- as.integer(snakemake@params[["min_hits"]])
PCT_GATE <- as.numeric(snakemake@params[["pct_gate"]])

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("reading eggNOG annotations")
ann <- read.delim(EGGNOG, comment.char = "#", header = FALSE, quote = "",
                  stringsAsFactors = FALSE)
colnames(ann)[c(1, 7, 13)] <- c("query", "COG_category", "KEGG_Pathway")
ann$protein <- sub(".*\\|", "", ann$query)

COG_NAMES <- c(
  J = "Translation, ribosomal structure and biogenesis",
  A = "RNA processing and modification",
  K = "Transcription",
  L = "Replication, recombination and repair",
  B = "Chromatin structure and dynamics",
  D = "Cell cycle control, cell division, chromosome partitioning",
  Y = "Nuclear structure",
  V = "Defense mechanisms",
  T = "Signal transduction mechanisms",
  M = "Cell wall/membrane/envelope biogenesis",
  N = "Cell motility",
  Z = "Cytoskeleton",
  W = "Extracellular structures",
  U = "Intracellular trafficking, secretion, and vesicular transport",
  O = "Posttranslational modification, protein turnover, chaperones",
  X = "Mobilome: prophages, transposons",
  C = "Energy production and conversion",
  G = "Carbohydrate transport and metabolism",
  E = "Amino acid transport and metabolism",
  F = "Nucleotide transport and metabolism",
  H = "Coenzyme transport and metabolism",
  I = "Lipid transport and metabolism",
  P = "Inorganic ion transport and metabolism",
  Q = "Secondary metabolites biosynthesis, transport and catabolism",
  R = "General function prediction only",
  S = "Function unknown")

kegg_names <- read.delim(KEGGNAME, stringsAsFactors = FALSE)
stopifnot(all(c("pathway_id", "name") %in% names(kegg_names)))
KEGG_NAMES <- setNames(kegg_names$name, sub("^[a-z]+", "", kegg_names$pathway_id))

cog_of <- local({
  v <- ann$COG_category
  v[v == "-" | is.na(v)] <- ""
  setNames(strsplit(v, ""), ann$protein)
})
kegg_cls <- read.delim(KEGGCLASS, comment.char = "#", stringsAsFactors = FALSE)
kegg_drop <- character(0)
if (DROP_NONBACT)  kegg_drop <- c(kegg_drop, kegg_cls$pathway_id[kegg_cls$is_non_bacterial == "yes"])
if (DROP_OVERVIEW) kegg_drop <- c(kegg_drop, kegg_cls$pathway_id[kegg_cls$is_overview == "yes"])
kegg_drop <- unique(sub("^ko", "map", kegg_drop))

kegg_known <- unique(sub("^ko", "map", kegg_cls$pathway_id))

kegg_of <- local({
  v <- ann$KEGG_Pathway
  v[v == "-" | is.na(v)] <- ""
  setNames(lapply(strsplit(v, ","), function(x) {
             t <- setdiff(x[startsWith(x, "map")], kegg_drop)
             if (DROP_RETIRED) intersect(t, kegg_known) else t
           }), ann$protein)
})
cat(sprintf("KEGG filter: %d maps dropped (%s)\n", length(kegg_drop),
            paste(c(if (DROP_NONBACT) "non-bacterial", if (DROP_OVERVIEW) "overview"),
                  collapse = " + ")))
cat(sprintf("KEGG terms retained on >=1 protein: %d\n",
            length(unique(unlist(kegg_of, use.names = FALSE)))))

CONTRASTS <- list(
  coculture_d2 = list(
    de     = "DE_coculture_d2.tsv",
    onoff  = "on_off_withC_d2_vs_alone_d2.tsv",
    up_grp = "withC_d2"))

rd <- function(f) read.delim(file.path(DE_DIR, f), stringsAsFactors = FALSE,
                             check.names = FALSE, quote = "")

enrich <- function(fg, universe, term_of, term_names, ontology) {
  fg <- intersect(fg, universe)
  annotated <- universe[lengths(term_of[universe]) > 0]
  fg <- intersect(fg, annotated)
  N <- length(annotated); n <- length(fg)
  if (n == 0) return(NULL)

  bg_terms <- unlist(term_of[annotated], use.names = FALSE)
  fg_terms <- unlist(term_of[fg],        use.names = FALSE)
  Kt <- table(bg_terms); kt <- table(fg_terms)

  terms <- names(Kt)[Kt >= MIN_TERM_SIZE]
  out <- lapply(terms, function(tm) {
    k <- if (tm %in% names(kt)) as.integer(kt[[tm]]) else 0L
    if (k < MIN_HITS) return(NULL)
    K <- as.integer(Kt[[tm]])
    m <- matrix(c(k, K - k, n - k, N - K - n + k), nrow = 2)
    ft <- fisher.test(m, alternative = "greater")
    hits <- fg[vapply(term_of[fg], function(x) tm %in% x, logical(1))]
    data.frame(ontology = ontology, term = tm,
               term_name = { key <- if (ontology == "KEGG") sub("^[a-z]+", "", tm) else tm
                             if (key %in% names(term_names)) term_names[[key]] else NA },
               k_fg = k, n_fg = n, K_bg = K, N_bg = N,
               fold_enrichment = round((k / n) / (K / N), 3),
               odds_ratio = round(unname(ft$estimate), 3),
               p_value = ft$p.value, proteins = paste(sort(hits), collapse = ","),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) return(NULL)
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out[order(out$p_value), ]
}

all_res <- list()
for (cname in names(CONTRASTS)) {
  cf <- CONTRASTS[[cname]]
  de <- rd(cf$de)

  oo <- if (is.na(cf$onoff)) NULL else rd(cf$onoff)
  universe <- unique(c(de$Protein, oo$Protein))

  de_up   <- de$Protein[de$sig_and_lfc & de$logFC > 0]
  de_down <- de$Protein[de$sig_and_lfc & de$logFC < 0]

  pct_sets <- function(thr) {
    if (is.null(oo)) return(list(up = character(0), down = character(0)))
    o <- oo[oo$pct_in_opposite_group > thr, ]
    list(up   = o$Protein[o$present_in  == cf$up_grp],
         down = o$Protein[o$absent_from == cf$up_grp])
  }
  pgate <- pct_sets(PCT_GATE)

  variants <- list(
    de_plus_onoff_p50  = list(up = union(de_up, pgate$up),
                              down = union(de_down, pgate$down)))

  for (vname in names(variants)) for (dir in c("up", "down")) {
    fg <- variants[[vname]][[dir]]
    for (ont in c("COG", "KEGG")) {
      r <- enrich(fg, universe,
                  if (ont == "COG") cog_of else kegg_of,
                  if (ont == "COG") COG_NAMES else KEGG_NAMES, ont)
      if (is.null(r)) next
      r <- cbind(contrast = cname, direction = dir, variant = vname, r)
      all_res[[length(all_res) + 1]] <- r
    }
    message(sprintf("  %-24s %-18s %-4s  fg=%d / universe=%d",
                    cname, vname, dir, length(intersect(fg, universe)), length(universe)))
  }
}
res <- do.call(rbind, all_res)

cols <- c("contrast", "direction", "variant", "ontology", "term", "term_name",
          "k_fg", "n_fg", "K_bg", "N_bg", "fold_enrichment", "odds_ratio",
          "p_value", "p_adj", "proteins")
res <- res[, cols]

for (ont in c("COG", "KEGG")) {
  f <- file.path(OUT_DIR, sprintf("%s_enrichment.tsv", ont))
  write.table(res[res$ontology == ont, ], f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", f)
}

SPLIT_OUT <- c(de_plus_onoff_p50 = PCT_GATE)
for (v in names(SPLIT_OUT)) {
  for (ont in c("COG", "KEGG")) {
    sub <- res[res$variant == v & res$ontology == ont, ]
    sub <- sub[order(sub$direction, sub$p_value), ]
    f <- file.path(OUT_DIR, sprintf("%s_enrichment_%s.tsv", ont, sub("^de_plus_onoff_", "onoff_", v)))
    con <- file(f, "w")
    writeLines(c(
      sprintf("# %s over-representation, contrast coculture_d2, variant '%s'", ont, v),
      "# foreground = DE-significant UNION ON/OFF calls. DE significance is the",
      "#   sig_and_lfc flag from 04_de_limma_MF6.R: adj.P < 0.05 AND |log2FC| >= 1.",
      "# ON/OFF calls are admitted only when their mean abundance exceeds the",
      sprintf("#   %gth percentile of the opposite group's detected abundances", SPLIT_OUT[[v]]),
      "#   (pct_in_opposite_group).",
      "# one-sided Fisher's exact, BH within this ontology x direction."), con)
    close(con)
    write.table(sub, f, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
    message("wrote ", f, "  (", nrow(sub), " terms, ",
            sum(sub$p_adj < 0.05, na.rm = TRUE), " significant)")
  }
}

message("\nsignificant terms (BH < 0.05):")
s <- res[res$p_adj < 0.05, ]
print(table(s$contrast, s$direction, s$ontology))
