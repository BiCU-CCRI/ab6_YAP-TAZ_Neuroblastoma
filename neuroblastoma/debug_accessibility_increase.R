#!/usr/bin/env Rscript

# Script to investigate CLB-Ma 24h accessibility increase
# Author: Claude Code Analysis

library(dplyr)
library(ggplot2)
library(corrplot)

cat("=== INVESTIGATING CLB-Ma 24h ACCESSIBILITY INCREASE ===\n\n")

# 1. Load the raw multiBigWigSummary data
cat("1. Loading raw ATAC-seq accessibility data...\n")
atac_scores <- read.table(
  "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript.tab",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Clean column names
colnames(atac_scores) <- gsub("'", "", colnames(atac_scores))
colnames(atac_scores) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores))

# Extract sample columns
sample_cols <- colnames(atac_scores)[4:ncol(atac_scores)]
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

print(sample_info)

# 2. Check library sizes (total accessibility scores per sample)
cat("\n2. Checking library sizes (total scores per sample)...\n")
accessibility_matrix <- as.matrix(atac_scores[, sample_cols])
library_sizes <- colSums(accessibility_matrix, na.rm = TRUE)

lib_size_df <- data.frame(
  sample = names(library_sizes),
  library_size = library_sizes
) %>%
  left_join(sample_info, by = "sample") %>%
  arrange(cell_line, timepoint, replicate)

print(lib_size_df)

# Calculate fold change in library size vs control
lib_size_summary <- lib_size_df %>%
  group_by(cell_line, timepoint) %>%
  summarise(mean_lib_size = mean(library_size), .groups = "drop") %>%
  arrange(cell_line, timepoint) %>%
  group_by(cell_line) %>%
  mutate(
    fc_vs_control = mean_lib_size / mean_lib_size[timepoint == "control"],
    log2fc_vs_control = log2(fc_vs_control)
  ) %>%
  ungroup()

cat("\nLibrary size fold changes vs control:\n")
print(lib_size_summary)

# 3. Check correlation between replicates
cat("\n3. Checking correlation between technical replicates...\n")
correlation_matrix <- cor(accessibility_matrix, use = "complete.obs")

# Extract replicate correlations
replicate_cors <- data.frame(
  condition = c("CLB_Ma_control", "CLB_Ma_24h", "CLB_Ma_48h", 
                "SK_N_SH_control", "SK_N_SH_24h", "SK_N_SH_48h"),
  correlation = c(
    correlation_matrix["CM_control_ATAC_S165169_REP1", "CM_control_ATAC_S165169_REP2"],
    correlation_matrix["CM_24h_ATAC_S165170_REP1", "CM_24h_ATAC_S165170_REP2"],
    correlation_matrix["CM_48h_ATAC_S165167_REP1", "CM_48h_ATAC_S165167_REP2"],
    correlation_matrix["SH_control_ATAC_S165166_REP1", "SH_control_ATAC_S165166_REP2"],
    correlation_matrix["SH_24h_ATAC_S165164_REP1", "SH_24h_ATAC_S165164_REP2"],
    correlation_matrix["SH_48h_ATAC_S165165_REP1", "SH_48h_ATAC_S165165_REP2"]
  )
)

print(replicate_cors)

# 4. Check distribution of accessibility scores
cat("\n4. Checking distribution of accessibility scores...\n")
# Calculate condition means
condition_means <- sample_info %>%
  group_by(cell_line, timepoint) %>%
  summarise(samples = list(sample), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    mean_accessibility = list(rowMeans(accessibility_matrix[, samples], na.rm = TRUE))
  ) %>%
  ungroup()

condition_matrix <- do.call(cbind, lapply(1:nrow(condition_means), function(i) {
  condition_means$mean_accessibility[[i]]
}))
colnames(condition_matrix) <- paste(condition_means$cell_line, condition_means$timepoint, sep = "_")

# Check median accessibility per condition
median_accessibility <- apply(condition_matrix, 2, median, na.rm = TRUE)
cat("Median accessibility per condition:\n")
print(median_accessibility)

# Calculate fold changes
clb_ma_24h_fc <- log2((condition_matrix[, "CLB_Ma_24h"] + 0.001) / 
                      (condition_matrix[, "CLB_Ma_control"] + 0.001))
clb_ma_48h_fc <- log2((condition_matrix[, "CLB_Ma_48h"] + 0.001) / 
                      (condition_matrix[, "CLB_Ma_control"] + 0.001))

cat(sprintf("\nCLB-Ma 24h: %.1f%% of regions show >0.5 log2FC increase\n",
            100 * sum(clb_ma_24h_fc > 0.5, na.rm = TRUE) / length(clb_ma_24h_fc)))
cat(sprintf("CLB-Ma 48h: %.1f%% of regions show >0.5 log2FC increase\n",
            100 * sum(clb_ma_48h_fc > 0.5, na.rm = TRUE) / length(clb_ma_48h_fc)))

# 5. Export results for further investigation
cat("\n5. Exporting diagnostic data...\n")
write.csv(lib_size_df, "~/workspace/neuroblastoma/temp_results/library_size_diagnostic.csv", row.names = FALSE)
write.csv(replicate_cors, "~/workspace/neuroblastoma/temp_results/replicate_correlations.csv", row.names = FALSE)

# Save accessibility summary
accessibility_summary <- data.frame(
  condition = colnames(condition_matrix),
  median_accessibility = median_accessibility,
  mean_accessibility = colMeans(condition_matrix, na.rm = TRUE),
  total_regions = nrow(condition_matrix)
)
write.csv(accessibility_summary, "~/workspace/neuroblastoma/temp_results/accessibility_summary.csv", row.names = FALSE)

cat("\n=== DIAGNOSTIC ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
cat(sprintf("- CLB-Ma 24h library size fold change: %.2f\n", 
            lib_size_summary$fc_vs_control[lib_size_summary$cell_line == "CLB_Ma" & lib_size_summary$timepoint == "24h"]))
cat(sprintf("- CLB-Ma 24h replicate correlation: %.3f\n", 
            replicate_cors$correlation[replicate_cors$condition == "CLB_Ma_24h"]))
cat("- Files saved to temp_results/ for further analysis\n")