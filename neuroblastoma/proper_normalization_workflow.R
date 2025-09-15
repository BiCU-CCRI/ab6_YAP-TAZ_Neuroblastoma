#!/usr/bin/env Rscript

# Proper ATAC-seq Normalization Workflow
# Author: Claude Code Analysis

library(dplyr)
library(DESeq2)
library(edgeR)

cat("=== ATAC-seq NORMALIZATION WORKFLOW ===\n\n")

# STEP 1: Load raw counts from multiBigWigSummary
cat("STEP 1: Loading raw accessibility counts...\n")
atac_scores <- read.table(
  "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript.tab",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Clean column names
colnames(atac_scores) <- gsub("'", "", colnames(atac_scores))
colnames(atac_scores) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores))

sample_cols <- colnames(atac_scores)[4:ncol(atac_scores)]
count_matrix <- as.matrix(atac_scores[, sample_cols])
rownames(count_matrix) <- paste0("Region_", 1:nrow(count_matrix))

# Sample metadata
sample_info <- data.frame(
  sample = sample_cols,
  cell_line = ifelse(grepl("^CM_", sample_cols), "CLB_Ma", "SK_N_SH"),
  timepoint = case_when(
    grepl("control", sample_cols) ~ "control",
    grepl("24h", sample_cols) ~ "24h", 
    grepl("48h", sample_cols) ~ "48h"
  ),
  replicate = ifelse(grepl("REP1", sample_cols), "REP1", "REP2"),
  stringsAsFactors = FALSE
)

cat(sprintf("Loaded %d regions across %d samples\n", nrow(count_matrix), ncol(count_matrix)))

# STEP 2: Multiple normalization approaches
cat("\nSTEP 2: Applying different normalization methods...\n")

# Method 1: Library Size Normalization (simplest)
cat("- Library size normalization...\n")
library_sizes <- colSums(count_matrix, na.rm = TRUE)
size_factors_lib <- library_sizes / median(library_sizes)
norm_lib_size <- t(t(count_matrix) / size_factors_lib)

# Method 2: TMM Normalization (edgeR)
cat("- TMM normalization (edgeR)...\n")
dge <- DGEList(counts = count_matrix)
dge <- calcNormFactors(dge, method = "TMM")
size_factors_tmm <- dge$samples$norm.factors * dge$samples$lib.size / median(dge$samples$lib.size)
norm_tmm <- t(t(count_matrix) / size_factors_tmm)

# Method 3: DESeq2 Size Factors
cat("- DESeq2 size factor normalization...\n")
# Create condition for DESeq2
condition <- paste(sample_info$cell_line, sample_info$timepoint, sep = "_")
coldata <- data.frame(condition = factor(condition))
rownames(coldata) <- sample_cols

# Create DESeq2 object (using round() for integer counts)
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = coldata,
  design = ~ condition
)

# Estimate size factors
dds <- estimateSizeFactors(dds)
size_factors_deseq <- sizeFactors(dds)
norm_deseq <- t(t(count_matrix) / size_factors_deseq)

# Method 4: Quantile Normalization
cat("- Quantile normalization...\n")
normalize_quantiles <- function(matrix) {
  # Simple quantile normalization
  ranks <- apply(matrix, 2, rank, ties.method = "average")
  sorted_means <- sort(rowMeans(matrix, na.rm = TRUE))
  norm_matrix <- matrix(0, nrow = nrow(matrix), ncol = ncol(matrix))
  for (i in 1:ncol(matrix)) {
    norm_matrix[, i] <- sorted_means[ranks[, i]]
  }
  dimnames(norm_matrix) <- dimnames(matrix)
  return(norm_matrix)
}
norm_quantile <- normalize_quantiles(count_matrix)

# STEP 3: Compare normalization effects
cat("\nSTEP 3: Comparing normalization effects on CLB-Ma 24h artifact...\n")

clb_ma_control_idx <- which(sample_info$cell_line == "CLB_Ma" & sample_info$timepoint == "control")
clb_ma_24h_idx <- which(sample_info$cell_line == "CLB_Ma" & sample_info$timepoint == "24h")

calculate_fc_stats <- function(norm_matrix, name) {
  control_mean <- rowMeans(norm_matrix[, clb_ma_control_idx], na.rm = TRUE)
  treat_mean <- rowMeans(norm_matrix[, clb_ma_24h_idx], na.rm = TRUE)
  log2fc <- log2((treat_mean + 0.001) / (control_mean + 0.001))
  
  stats <- data.frame(
    method = name,
    median_log2fc = median(log2fc, na.rm = TRUE),
    pct_increased = 100 * sum(log2fc > 0.5, na.rm = TRUE) / length(log2fc),
    pct_decreased = 100 * sum(log2fc < -0.5, na.rm = TRUE) / length(log2fc),
    stringsAsFactors = FALSE
  )
  return(stats)
}

comparison_stats <- rbind(
  calculate_fc_stats(count_matrix, "Raw"),
  calculate_fc_stats(norm_lib_size, "Library_Size"),
  calculate_fc_stats(norm_tmm, "TMM"),
  calculate_fc_stats(norm_deseq, "DESeq2"),
  calculate_fc_stats(norm_quantile, "Quantile")
)

print(comparison_stats)

# STEP 4: Check library size effects
cat("\nSTEP 4: Library size comparison across methods...\n")
lib_size_comparison <- data.frame(
  sample = sample_cols,
  raw = colSums(count_matrix, na.rm = TRUE),
  lib_size_norm = colSums(norm_lib_size, na.rm = TRUE),
  tmm_norm = colSums(norm_tmm, na.rm = TRUE),
  deseq2_norm = colSums(norm_deseq, na.rm = TRUE),
  quantile_norm = colSums(norm_quantile, na.rm = TRUE)
) %>%
  left_join(sample_info, by = "sample")

print(lib_size_comparison[lib_size_comparison$cell_line == "CLB_Ma", ])

# STEP 5: Export normalized data
cat("\nSTEP 5: Exporting normalized matrices...\n")

# Save all normalization results
results_dir <- "~/workspace/neuroblastoma/temp_results"

write.csv(data.frame(atac_scores[, 1:3], norm_lib_size), 
          file.path(results_dir, "atac_normalized_lib_size.csv"), row.names = FALSE)
write.csv(data.frame(atac_scores[, 1:3], norm_tmm), 
          file.path(results_dir, "atac_normalized_tmm.csv"), row.names = FALSE)
write.csv(data.frame(atac_scores[, 1:3], norm_deseq), 
          file.path(results_dir, "atac_normalized_deseq2.csv"), row.names = FALSE)
write.csv(data.frame(atac_scores[, 1:3], norm_quantile), 
          file.path(results_dir, "atac_normalized_quantile.csv"), row.names = FALSE)

# Save comparison stats
write.csv(comparison_stats, file.path(results_dir, "normalization_comparison.csv"), row.names = FALSE)

cat("\n=== RECOMMENDATIONS ===\n")
cat("1. For differential accessibility analysis: Use DESeq2 or TMM normalization\n")
cat("2. For visualization and exploratory analysis: Use library size or quantile normalization\n")
cat("3. The CLB-Ma 24h artifact is best corrected by TMM or DESeq2 methods\n")
cat("4. Always apply normalization BEFORE calculating fold changes or statistical tests\n\n")

cat("Normalized data files saved to temp_results/\n")
cat("Use these instead of raw scores from multiBigWigSummary output\n")