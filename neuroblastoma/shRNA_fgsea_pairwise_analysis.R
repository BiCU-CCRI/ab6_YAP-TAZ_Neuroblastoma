# shRNA-seq Pairwise fGSEA Analysis
# Analyzing DOX vs CTRL pairs without replicates
# Author: Claude Code
# Date: 2025-09-15

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(hypeR)
  library(dplyr)
  library(ggplot2)
  library(here)
  library(openxlsx2)
})

# Configuration
PATHS <- list(
  data_dir = here::here("neuroblastoma", "data", "RNAseq"),
  results_dir = here::here("neuroblastoma", "results", "RNA-seq_shRNA_pairwise"),
  resources_dir = here::here("neuroblastoma", "resources")
)

# Create results directory
if(!dir.exists(PATHS$results_dir)) dir.create(PATHS$results_dir, recursive = TRUE)

# Load utility functions
source(here::here("neuroblastoma/resources/UtilityScriptsRNA-seq.R"))

# Load annotation data
annotation_data <- read.table(
  file = file.path(PATHS$data_dir, "rnaseq_deseq_global_annotation_gene.tsv"),
  sep = "\t",
  header = TRUE
)
colnames(annotation_data)[c(5, 7)] <- c("ensembl_id", "gene_symbol")

# Load shRNA DESeq2 object
message("Loading shRNA-seq data...")
load(file.path(PATHS$data_dir, "deseq2_shRNA.dds.RData"))

# Examine the data structure
message("Examining data structure...")
print("Sample metadata:")
print(colData(dds))
print("Design formula:")
print(design(dds))
print("Sample names:")
print(colnames(dds))

# Identify sample pairs based on metadata
# Assuming Group1 contains the treatment types (LUC, TGFP, YAP)
# and samples alternate between CTRL and DOX
sample_info <- as.data.frame(colData(dds))
sample_info$sample_name <- rownames(sample_info)

# Print sample information to understand structure
print("Full sample information:")
print(sample_info)

# Basic filtering and normalization
message("Filtering and normalizing data...")

# Filter lowly expressed genes (keep genes with ≥10 counts in ≥2 samples)
keep <- rowSums(counts(dds) >= 10) >= 2
dds_filtered <- dds[keep,]

# Estimate size factors for normalization
dds_filtered <- estimateSizeFactors(dds_filtered)

# Get normalized counts
normalized_counts <- counts(dds_filtered, normalized = TRUE)

# Add gene symbols as rownames for easier interpretation
gene_symbols <- annotation_data$gene_symbol[match(rownames(normalized_counts), annotation_data$ensembl_id)]
# Handle missing symbols
gene_symbols[is.na(gene_symbols)] <- rownames(normalized_counts)[is.na(gene_symbols)]

# Handle duplicate gene symbols by keeping the gene with highest mean expression
mean_expr <- rowMeans(normalized_counts)
df_temp <- data.frame(
  ensembl_id = rownames(normalized_counts),
  gene_symbol = gene_symbols, 
  mean_expr = mean_expr,
  stringsAsFactors = FALSE
)

# Remove duplicates by keeping the gene with highest mean expression
df_temp <- df_temp %>%
  group_by(gene_symbol) %>%
  slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup()

# Filter normalized counts to remove duplicates and set gene symbols as rownames
normalized_counts <- normalized_counts[df_temp$ensembl_id, ]
rownames(normalized_counts) <- df_temp$gene_symbol

# Load gene sets for fGSEA
message("Loading gene sets...")
gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean = TRUE)
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean = TRUE)
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean = TRUE)
gs_C6_onco <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C6"), clean = TRUE)

# Custom gene sets
gs_wang_hippo <- list(wang_hippo_set = c("CCN1", "CCN2", "AMOTL2", "ANKRD1", "IGFBP3", "F3", "FJX1", "NUAK2", "LATS2", "CRIM1", "GADD45A",
                                         "TGFB2", "PTPN14", "NT5E", "FOXF2", "AXL", "DOCK5", "ASAP1", "RBMS3", "MYOF", "ARHGEF17", "CCDC80"))

# Combine gene sets for analysis
gene_set_list <- list(
  Hallmark_EMT = gs_hallmark[["genesets"]][["Epithelial Mesenchymal Transition"]],
  Hallmark_Hippo = gs_hallmark[["genesets"]][["Hippo Signaling"]],
  GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
  REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
  C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
  Hippo_Wang = gs_wang_hippo$wang_hippo_set
)

# Function to calculate log2FC and run fGSEA for a pair of samples
analyze_sample_pair <- function(ctrl_sample, dox_sample, pair_name, normalized_counts, gene_sets) {
  message(paste("Analyzing pair:", pair_name))
  
  # Get normalized counts for the pair
  ctrl_counts <- normalized_counts[, ctrl_sample]
  dox_counts <- normalized_counts[, dox_sample]
  
  # Calculate log2FC with pseudocount to avoid division by zero
  pseudocount <- 1
  log2fc <- log2((dox_counts + pseudocount) / (ctrl_counts + pseudocount))
  
  # Remove genes with NA or infinite values
  valid_genes <- is.finite(log2fc) & !is.na(log2fc)
  log2fc <- log2fc[valid_genes]
  
  # Add small random jitter to break ties
  set.seed(42)  # for reproducibility
  log2fc <- log2fc + runif(length(log2fc), -1e-8, 1e-8)
  
  # Sort genes by log2FC (descending order for fGSEA)
  log2fc_sorted <- sort(log2fc, decreasing = TRUE)
  
  # Run fGSEA
  fgsea_results <- fgsea(
    pathways = gene_sets,
    stats = log2fc_sorted,
    minSize = 15,
    maxSize = 2500
  )
  
  # Add pair information to results
  fgsea_results$pair <- pair_name
  fgsea_results$ctrl_sample <- ctrl_sample
  fgsea_results$dox_sample <- dox_sample
  
  return(list(
    fgsea_results = fgsea_results,
    log2fc_ranked = log2fc_sorted,
    pair_name = pair_name
  ))
}

# Function to create enrichment plots
create_enrichment_plots <- function(analysis_result, gene_sets, output_dir) {
  fgsea_res <- analysis_result$fgsea_results
  log2fc_ranked <- analysis_result$log2fc_ranked
  pair_name <- analysis_result$pair_name
  
  # Create plots for significant pathways
  significant_pathways <- fgsea_res[fgsea_res$padj < 0.05, ]
  
  if(nrow(significant_pathways) > 0) {
    for(i in 1:nrow(significant_pathways)) {
      pathway_name <- significant_pathways$pathway[i]
      pathway_genes <- gene_sets[[pathway_name]]
      
      if(!is.null(pathway_genes) && length(pathway_genes) > 0) {
        tryCatch({
          # Create enrichment plot
          p <- plotEnrichment(pathway_genes, log2fc_ranked) +
            labs(
              title = paste(pair_name, "-", pathway_name),
              subtitle = paste("NES =", round(significant_pathways$NES[i], 3),
                             "| padj =", format(significant_pathways$padj[i], scientific = TRUE, digits = 3))
            ) +
            theme_minimal()
          
          # Save plot
          filename <- file.path(output_dir, paste0(pair_name, "_", gsub("[^A-Za-z0-9]", "_", pathway_name), "_enrichment.pdf"))
          ggsave(filename, p, width = 10, height = 8)
        }, error = function(e) {
          message(paste("Warning: Could not create plot for", pathway_name, ":", e$message))
        })
      }
    }
  } else {
    message(paste("No significant pathways found for", pair_name))
  }
}

# Main analysis workflow
message("Starting pairwise analysis...")

# Based on typical shRNA experimental design, assuming we have:
# Samples 1,3,5 = CTRL and Samples 2,4,6 = DOX for LUC, TGFP, YAP respectively
# This needs to be verified based on actual sample metadata

# Let's first check the actual sample organization
sample_pairs <- list()
ctrl_samples <- c()
dox_samples <- c()

# Try to identify pairs based on sample names and metadata
# This is a placeholder - adjust based on actual data structure
if("Group1" %in% colnames(sample_info)) {
  unique_groups <- unique(sample_info$Group1)
  
  for(group in unique_groups) {
    group_samples <- sample_info[sample_info$Group1 == group, ]
    
    # Assuming we have one CTRL and one DOX per group
    # Adjust this logic based on actual naming convention
    if(nrow(group_samples) == 2) {
      # Try to identify which is CTRL vs DOX based on sample names or other metadata
      sample_names <- rownames(group_samples)
      
      # Simple heuristic - adjust as needed
      ctrl_idx <- grep("ctrl|control|CTRL|CONTROL", sample_names, ignore.case = TRUE)
      dox_idx <- grep("dox|DOX|treat|TREAT", sample_names, ignore.case = TRUE)
      
      if(length(ctrl_idx) == 0) ctrl_idx <- 1  # fallback
      if(length(dox_idx) == 0) dox_idx <- 2    # fallback
      
      sample_pairs[[paste0(group, "_DOX_vs_CTRL")]] <- list(
        ctrl = sample_names[ctrl_idx[1]],
        dox = sample_names[dox_idx[1]]
      )
    }
  }
} else {
  # Fallback: assume alternating pattern
  message("Using fallback sample pairing based on position...")
  sample_names <- colnames(dds_filtered)
  
  # Assuming 6 samples: 1,2 = LUC pair; 3,4 = TGFP pair; 5,6 = YAP pair
  if(length(sample_names) >= 6) {
    sample_pairs <- list(
      "LUC_DOX_vs_CTRL" = list(ctrl = sample_names[1], dox = sample_names[2]),
      "TGFP_DOX_vs_CTRL" = list(ctrl = sample_names[3], dox = sample_names[4]),
      "YAP_DOX_vs_CTRL" = list(ctrl = sample_names[5], dox = sample_names[6])
    )
  }
}

# Display identified pairs for verification
message("Identified sample pairs:")
print(sample_pairs)

# Run analysis for each pair
all_results <- list()
combined_fgsea <- data.frame()

for(pair_name in names(sample_pairs)) {
  pair_info <- sample_pairs[[pair_name]]
  
  # Run analysis
  result <- analyze_sample_pair(
    ctrl_sample = pair_info$ctrl,
    dox_sample = pair_info$dox,
    pair_name = pair_name,
    normalized_counts = normalized_counts,
    gene_sets = gene_set_list
  )
  
  all_results[[pair_name]] <- result
  
  # Combine results
  if(nrow(combined_fgsea) == 0) {
    combined_fgsea <- result$fgsea_results
  } else {
    combined_fgsea <- rbind(combined_fgsea, result$fgsea_results)
  }
  
  # Create enrichment plots
  create_enrichment_plots(result, gene_set_list, PATHS$results_dir)
  
  # Save individual results (convert list columns to character for Excel compatibility)
  fgsea_for_excel <- result$fgsea_results
  fgsea_for_excel$leadingEdge <- sapply(fgsea_for_excel$leadingEdge, function(x) paste(x, collapse = "; "))
  
  openxlsx2::write_xlsx(
    fgsea_for_excel,
    file.path(PATHS$results_dir, paste0(pair_name, "_fgsea_results.xlsx"))
  )
}

# Save combined results (convert list columns to character for Excel compatibility)
combined_fgsea_for_excel <- combined_fgsea
if(nrow(combined_fgsea_for_excel) > 0) {
  combined_fgsea_for_excel$leadingEdge <- sapply(combined_fgsea_for_excel$leadingEdge, function(x) paste(x, collapse = "; "))
}
openxlsx2::write_xlsx(combined_fgsea_for_excel, file.path(PATHS$results_dir, "combined_fgsea_results.xlsx"))

# Create summary visualization
message("Creating summary visualization...")

# Summary heatmap of NES values
library(reshape2)
library(pheatmap)

# Prepare data for heatmap
heatmap_data <- combined_fgsea %>%
  select(pathway, pair, NES, padj) %>%
  filter(padj < 0.05) %>%  # Only significant pathways
  dcast(pathway ~ pair, value.var = "NES", fill = 0)

if(nrow(heatmap_data) > 1 && ncol(heatmap_data) > 2) {  # Need at least 2 pathways and 2 pairs
  rownames(heatmap_data) <- heatmap_data$pathway
  heatmap_data <- heatmap_data[, -1, drop = FALSE]
  
  # Determine clustering options based on data dimensions
  cluster_rows <- nrow(heatmap_data) > 1
  cluster_cols <- ncol(heatmap_data) > 1
  
  # Create heatmap
  pheatmap(
    heatmap_data,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "NES Heatmap - Significant Pathways (padj < 0.05)",
    filename = file.path(PATHS$results_dir, "fgsea_NES_heatmap.pdf"),
    width = 10,
    height = 8
  )
  message("Heatmap created successfully")
} else {
  message("Insufficient data for heatmap generation (need at least 2 significant pathways)")
}

# Create summary table
summary_table <- combined_fgsea %>%
  filter(padj < 0.05) %>%
  arrange(pair, padj) %>%
  select(pair, pathway, NES, pval, padj, size)

openxlsx2::write_xlsx(summary_table, file.path(PATHS$results_dir, "significant_pathways_summary.xlsx"))

message("Analysis complete!")
message(paste("Results saved to:", PATHS$results_dir))
message(paste("Number of sample pairs analyzed:", length(sample_pairs)))
message(paste("Total significant pathways (padj < 0.05):", nrow(summary_table)))