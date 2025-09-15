#!/usr/bin/env Rscript

# Visualization and correction strategies for CLB-Ma 24h accessibility increase
# Author: Claude Code Analysis

library(dplyr)
library(ggplot2)
library(tidyr)

# Load diagnostic data
lib_size_df <- read.csv("~/workspace/neuroblastoma/temp_results/library_size_diagnostic.csv")
accessibility_summary <- read.csv("~/workspace/neuroblastoma/temp_results/accessibility_summary.csv")

# Create visualization directory
results_dir <- "~/workspace/neuroblastoma/results/diagnostic_plots"
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 1. Library size visualization
p1 <- ggplot(lib_size_df, aes(x = timepoint, y = library_size, fill = cell_line)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.1f", library_size)), 
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  facet_wrap(~replicate) +
  labs(
    title = "Library Sizes Across Conditions",
    subtitle = "CLB-Ma 24h shows ~13% increase vs control",
    x = "Treatment Timepoint",
    y = "Total Accessibility Score",
    fill = "Cell Line"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "library_size_comparison.pdf"), p1, width = 10, height = 6)
print(p1)

# 2. Fold change visualization
lib_size_summary <- lib_size_df %>%
  group_by(cell_line, timepoint) %>%
  summarise(mean_lib_size = mean(library_size), .groups = "drop") %>%
  arrange(cell_line, timepoint) %>%
  group_by(cell_line) %>%
  mutate(
    fc_vs_control = mean_lib_size / mean_lib_size[timepoint == "control"],
    log2fc_vs_control = log2(fc_vs_control)
  ) %>%
  ungroup() %>%
  filter(timepoint != "control")

p2 <- ggplot(lib_size_summary, aes(x = timepoint, y = log2fc_vs_control, fill = cell_line)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.2f", log2fc_vs_control)), 
            position = position_dodge(width = 0.9), vjust = -0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "Library Size Changes vs Control",
    subtitle = "CLB-Ma 24h shows significant increase (+0.18 log2FC)",
    x = "Treatment Timepoint",
    y = "Log2 Fold Change vs Control",
    fill = "Cell Line"
  ) +
  theme_minimal()

ggsave(file.path(results_dir, "library_size_fold_changes.pdf"), p2, width = 8, height = 6)
print(p2)

# 3. Propose correction strategy
cat("\n=== CORRECTION STRATEGY ===\n")
cat("The CLB-Ma 24h samples show a systematic ~13% increase in total accessibility.\n")
cat("This affects both replicates consistently, suggesting a technical artifact.\n\n")

cat("Recommended correction approaches:\n")
cat("1. QUANTILE NORMALIZATION: Normalize samples to same distribution\n")
cat("2. TMM NORMALIZATION: Use trimmed mean of M-values (edgeR/DESeq2)\n")
cat("3. LIBRARY SIZE NORMALIZATION: Scale by total counts\n")
cat("4. SPIKE-IN NORMALIZATION: If spike-ins were used during sample prep\n\n")

# 4. Apply library size normalization as example
cat("Applying library size normalization as example...\n")

# Load raw data for normalization
atac_scores <- read.table(
  "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript.tab",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# Clean column names
colnames(atac_scores) <- gsub("'", "", colnames(atac_scores))
colnames(atac_scores) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores))

sample_cols <- colnames(atac_scores)[4:ncol(atac_scores)]
accessibility_matrix <- as.matrix(atac_scores[, sample_cols])

# Calculate normalization factors
library_sizes <- colSums(accessibility_matrix, na.rm = TRUE)
size_factors <- library_sizes / median(library_sizes)

# Apply normalization
normalized_matrix <- t(t(accessibility_matrix) / size_factors)
colnames(normalized_matrix) <- sample_cols

# Check effect of normalization
norm_lib_sizes <- colSums(normalized_matrix, na.rm = TRUE)
cat("Library sizes after normalization:\n")
norm_summary <- data.frame(
  sample = names(norm_lib_sizes),
  original = library_sizes,
  normalized = norm_lib_sizes,
  size_factor = size_factors
) %>%
  left_join(lib_size_df[, c("sample", "cell_line", "timepoint")], by = "sample")

print(norm_summary)

# 5. Compare original vs normalized results for CLB-Ma 24h
cat("\n=== COMPARISON: ORIGINAL vs NORMALIZED ===\n")

# Calculate condition means for normalized data
sample_info <- data.frame(
  sample = sample_cols,
  cell_line = ifelse(grepl("^CM_", sample_cols), "CLB_Ma", "SK_N_SH"),
  timepoint = case_when(
    grepl("control", sample_cols) ~ "control",
    grepl("24h", sample_cols) ~ "24h", 
    grepl("48h", sample_cols) ~ "48h"
  ),
  stringsAsFactors = FALSE
)

clb_ma_control_samples <- sample_info$sample[sample_info$cell_line == "CLB_Ma" & sample_info$timepoint == "control"]
clb_ma_24h_samples <- sample_info$sample[sample_info$cell_line == "CLB_Ma" & sample_info$timepoint == "24h"]

# Original fold changes
original_control_mean <- rowMeans(accessibility_matrix[, clb_ma_control_samples], na.rm = TRUE)
original_24h_mean <- rowMeans(accessibility_matrix[, clb_ma_24h_samples], na.rm = TRUE)
original_fc <- log2((original_24h_mean + 0.001) / (original_control_mean + 0.001))

# Normalized fold changes  
norm_control_mean <- rowMeans(normalized_matrix[, clb_ma_control_samples], na.rm = TRUE)
norm_24h_mean <- rowMeans(normalized_matrix[, clb_ma_24h_samples], na.rm = TRUE)
norm_fc <- log2((norm_24h_mean + 0.001) / (norm_control_mean + 0.001))

comparison_summary <- data.frame(
  metric = c("Median log2FC", "% regions >0.5 log2FC", "% regions <-0.5 log2FC"),
  original = c(
    median(original_fc, na.rm = TRUE),
    100 * sum(original_fc > 0.5, na.rm = TRUE) / length(original_fc),
    100 * sum(original_fc < -0.5, na.rm = TRUE) / length(original_fc)
  ),
  normalized = c(
    median(norm_fc, na.rm = TRUE),
    100 * sum(norm_fc > 0.5, na.rm = TRUE) / length(norm_fc),
    100 * sum(norm_fc < -0.5, na.rm = TRUE) / length(norm_fc)
  )
)

print(comparison_summary)

# Save normalized data
write.csv(data.frame(atac_scores[, 1:3], normalized_matrix), 
          "~/workspace/neuroblastoma/temp_results/normalized_accessibility_scores.csv", 
          row.names = FALSE)

cat("\nNormalized data saved to temp_results/normalized_accessibility_scores.csv\n")
cat("Diagnostic plots saved to results/diagnostic_plots/\n")