log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(topGO)
library(dplyr)
library(tidyr)

# Snakemake in/out and params
dds_file_path   <- snakemake@input[["deseq2_rds"]]
all_genes_path  <- snakemake@input[["all_annotations"]]
up_clust_path <- if (!is.null(snakemake@input[["up_clusters"]])) {
  snakemake@input[["up_clusters"]]
} else {
  NULL
}

down_clust_path <- if (!is.null(snakemake@input[["down_clusters"]])) {
  snakemake@input[["down_clusters"]]
} else {
  NULL
}

single_file_path <- if (!is.null(snakemake@input[["single_file"]])) {
  snakemake@input[["single_file"]]
} else {
  NULL
}
print(up_clust_path)
print(down_clust_path)
print(single_file_path)
out_tsv <- snakemake@output[["enrichment_table"]]

min_membership <- as.numeric(snakemake@params[["minMembership"]])
min_node_size <- as.numeric(snakemake@params[["minNodeSize"]])
min_gene_count <- as.numeric(snakemake@params[["minGeneCount"]])

# Load and filter
all_genes_df <- read.csv(all_genes_path, stringsAsFactors = FALSE, sep = "\t")
dds <- readRDS(dds_file_path)
expressed_genes <- rownames(dds)

# Filter Universe
universe_df <- all_genes_df %>% filter(Gene %in% expressed_genes)
universe_vector <- unique(universe_df$Gene)

# GENE-TO-GO mapping
message("Building Gene-to-GO mapping...")
gene_go_long <- universe_df %>%
  select(Gene, GO) %>%
  mutate(GO = trimws(GO)) %>%
  filter(!is.na(GO) & GO != "") %>%
  separate_rows(GO, sep = " +") %>%
  filter(GO != "")

gene2go <- split(gene_go_long$GO, gene_go_long$Gene)

# Prepare clusters
if (!is.null(up_clust_path) &&
    !is.null(down_clust_path)) {

  up_df <- read.csv(up_clust_path,
                    stringsAsFactors = FALSE, sep = "\t") %>%
    mutate(Cluster_Label = paste0("Up_C", Cluster))

  down_df <- read.csv(down_clust_path,
                      stringsAsFactors = FALSE, sep = "\t") %>%
    mutate(Cluster_Label = paste0("Down_C", Cluster))

  merged_clusters <- bind_rows(up_df, down_df) %>%
    filter(Membership_Score >= min_membership)

} else if (!is.null(single_file_path)) {

  merged_clusters <- read.csv(single_file_path,
                              stringsAsFactors = FALSE, sep = "\t") %>%
    mutate(Cluster_Label = paste0("", Cluster))

} else {
  stop("Need either both up_clust_path and down_clust_path, or single_file_path.")
}

unique_clusters <- unique(merged_clusters$Cluster_Label)

# TopGO func
run_topgo_namespace <- function(ontology_type, gene_factor, gene2go_map) {

  go_data <- new("topGOdata",
                 description = paste("Analysis for", ontology_type),
                 ontology = ontology_type,
                 allGenes = gene_factor,
                 nodeSize = min_node_size,
                 annot = annFUN.gene2GO,
                 gene2GO = gene2go_map)

  result_weight01 <- runTest(go_data,
                             algorithm = "weight01",
                             statistic = "fisher")

  # Get Top 20 terms
  res_table <- GenTable(go_data,
                        weight01Fisher = result_weight01,
                        orderBy = "weight01Fisher",
                        ranksOf = "weight01Fisher",
                        topNodes = 20,
                        numChar = 1000)
  # Calculate -log10
  # Convert scientific text to numeric form
  p_clean <- gsub("[< ]", "", res_table$weight01Fisher)
  # Convert to numeric
  p_numeric <- as.numeric(p_clean)
  neg_log10_p <- -log10(p_numeric)
  neg_log10_p[is.infinite(neg_log10_p)] <- 300  # Cap infinity

  res_table$neg_log10_p <- neg_log10_p

  sig_genes <- names(gene_factor)[as.integer(as.character(gene_factor)) == 1]

  # For each GO term in the results, extract the genes
  res_table$geneID <- sapply(res_table$GO.ID, function(go_id) {

    genes_in_term <- genesInTerm(go_data, go_id)[[1]]
    genes_overlap <- intersect(genes_in_term, sig_genes)
    paste(genes_overlap, collapse = "|")
  })

  return(res_table)
}

all_results_list <- list()

for (cluster_name in unique_clusters) {
  message(paste("Processing Cluster:", cluster_name))

  # Define target genes
  target_genes <- merged_clusters %>%
    filter(Cluster_Label == cluster_name) %>%
    pull(Gene)

  gene_factor <- factor(as.integer(universe_vector %in% target_genes),
                        levels = c(0, 1))
  names(gene_factor) <- universe_vector

  # Skip if cluster is too small
  if (sum(as.numeric(as.character(gene_factor))) < min_gene_count) {
    message(paste("Skipping", cluster_name, "- Too few genes."))
    next
  }

  # Loop Ontologies
  for (ontology in c("BP", "CC", "MF")) {
    tryCatch({
      # Run Analysis
      res_tbl <- run_topgo_namespace(ontology, gene_factor, gene2go)
      # Add Metadata columns
      res_tbl$Cluster <- cluster_name
      res_tbl$Ontology <- ontology
      # Append to master list
      all_results_list[[paste(cluster_name, ontology, sep = "_")]] <- res_tbl
    }, error = function(e) {
      message(paste("Error in", cluster_name, ontology, ":", e$message))
    })
  }
}

# Save table
final_master_table <- bind_rows(all_results_list)
final_master_table <- final_master_table %>%
  arrange(desc(neg_log10_p)) %>%
  select(Cluster, Ontology, GO.ID, Term,
         weight01Fisher, Significant, Annotated, everything())

write.table(final_master_table, file = out_tsv,
            row.names = FALSE, sep = "\t", quote = FALSE)