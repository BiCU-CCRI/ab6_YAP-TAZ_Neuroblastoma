library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggpubr)

library(GenomicRanges)

# Setup paths
data_dir <- "~/workspace/neuroblastoma/data"
results_dir <- "~/workspace/neuroblastoma/results/ATAC-seq_drug_treatment/unique_enhancers"

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


# Define file paths and metadata
files_to_process <- list(
  SKN_A_u = list(
    path = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_skn_ADRu.tab",
    cell_line = "SK_N_SH",
    phenotype = "ADR",
    description = "SK-N-SH ADR-unique super enhancers"
  ),
  SKN_M_u = list(
    path = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_skn_MESu.tab",
    cell_line = "SK_N_SH", 
    phenotype = "MES",
    description = "SK-N-SH MES-unique super enhancers"
  ),
  CLMB_A_u = list(
    path = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_clbm_ADRu.tab",
    cell_line = "CLB_Ma",
    phenotype = "ADR", 
    description = "CLB-Ma ADR-unique super enhancers"
  ),
  CLMB_M_u = list(
    path = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_clbm_MESu.tab",
    cell_line = "CLB_Ma",
    phenotype = "MES",
    description = "CLB-Ma MES-unique super enhancers"
  )
)

cat("Loading raw ATAC-seq accessibility scores for unique super enhancers...\n")
list_of_SEs <- list()

# Load and process each file
for (file_name in names(files_to_process)) {
  file_info <- files_to_process[[file_name]]
  cat(sprintf("Loading %s...\n", file_info$description))
  
  atac_scores_raw <- read.table(
    file_info$path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  # Clean column names
  colnames(atac_scores_raw) <- gsub("'", "", colnames(atac_scores_raw))
  colnames(atac_scores_raw) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores_raw))
  
  # Add metadata
  atac_scores_raw$file_type <- file_name
  atac_scores_raw$cell_line <- file_info$cell_line
  atac_scores_raw$phenotype <- file_info$phenotype
  
  list_of_SEs[[file_name]] <- atac_scores_raw
  
  cat(sprintf("Loaded %d unique %s super enhancers for %s\n", 
              nrow(atac_scores_raw), file_info$phenotype, file_info$cell_line))
}

# Summary of loaded data
total_regions <- sum(sapply(list_of_SEs, nrow))
cat(sprintf("\nTotal unique super enhancer regions loaded: %d\n", total_regions))
cat("Breakdown by type:\n")
for (file_name in names(list_of_SEs)) {
  cat(sprintf("  %s: %d regions\n", files_to_process[[file_name]]$description, nrow(list_of_SEs[[file_name]])))
}








cat("Applying genome-wide normalization to correct for sequencing depth bias...\n")

# Load genome-wide size factors calculated from chromosome 1 sampling
size_factors_file <- "~/workspace/neuroblastoma/temp_results/genome_wide_size_factors.csv"

if (file.exists(size_factors_file)) {
  # Load pre-calculated size factors
  sample_size_factors <- read.csv(size_factors_file)
  cat("Using pre-calculated genome-wide size factors\n")
  
} else {
  # Calculate genome-wide size factors on-the-fly using rtracklayer
  cat("Calculating genome-wide size factors from bigWig files...\n")
  
  library(rtracklayer)
  
  # Get bigWig files
  bigwig_dir <- file.path(data_dir, "ATACseq/bigWigs")
  bigwig_files <- list.files(bigwig_dir, pattern = "\\.bigWig$", full.names = TRUE)
  
  # Create sampling bins on chromosome 1
  chr1_bins <- GRanges(
    seqnames = "1",
    ranges = IRanges(
      start = seq(1, 250000000, by = 100000),
      width = 100000
    )
  )
  
  # Calculate coverage for each sample
  sample_coverage <- data.frame(
    sample = gsub("\\.mLb\\.clN\\.bigWig", "", basename(bigwig_files)),
    total_coverage = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(bigwig_files)) {
    coverage_data <- import(bigwig_files[i], format = "BigWig", which = chr1_bins)
    sample_coverage$total_coverage[i] <- sum(coverage_data$score * width(coverage_data), na.rm = TRUE)
  }
  
  # Calculate size factors
  median_coverage <- median(sample_coverage$total_coverage, na.rm = TRUE)
  sample_size_factors <- sample_coverage %>%
    mutate(
      cell_line = ifelse(grepl("^CM_", sample), "CLB_Ma", "SK_N_SH"),
      timepoint = case_when(
        grepl("control", sample) ~ "control",
        grepl("24h", sample) ~ "24h", 
        grepl("48h", sample) ~ "48h"
      ),
      replicate = ifelse(grepl("REP1", sample), "REP1", "REP2"),
      size_factor = total_coverage / median_coverage
    )
  
  # Save for future use
  write.csv(sample_size_factors, size_factors_file, row.names = FALSE)
}

# Display size factor summary
cat("\nGenome-wide size factors:\n")
size_factor_summary <- sample_size_factors %>%
  group_by(cell_line, timepoint) %>%
  summarise(mean_size_factor = mean(size_factor, na.rm = TRUE), .groups = "drop")
print(size_factor_summary)

# Apply genome-wide normalization to each dataset
normalized_list_of_SEs <- list()

for (file_name in names(list_of_SEs)) {
  cat(sprintf("\nApplying normalization to %s...\n", files_to_process[[file_name]]$description))
  
  atac_scores_raw <- list_of_SEs[[file_name]]
  
  # Get sample columns (excluding coordinate and metadata columns)
  sample_cols <- colnames(atac_scores_raw)[4:(ncol(atac_scores_raw)-3)]  # Exclude chr, start, end, file_type, cell_line, phenotype
  atac_matrix_raw <- as.matrix(atac_scores_raw[, sample_cols])
  
  # Create normalized matrix
  atac_matrix_normalized <- atac_matrix_raw
  for (i in seq_along(sample_cols)) {
    sample_name <- sample_cols[i]
    size_factor <- sample_size_factors$size_factor[sample_size_factors$sample == sample_name]
    
    if (length(size_factor) == 1 && !is.na(size_factor)) {
      atac_matrix_normalized[, i] <- atac_matrix_raw[, i] / size_factor
    }
  }
  
  # Create final normalized data frame
  atac_scores_normalized <- data.frame(
    atac_scores_raw[, 1:3],  # Keep coordinate columns
    atac_matrix_normalized,
    file_type = atac_scores_raw$file_type,
    cell_line = atac_scores_raw$cell_line,
    phenotype = atac_scores_raw$phenotype,
    stringsAsFactors = FALSE
  )
  
  # Validation: Calculate normalization effect for CLB-Ma samples if present
  if (files_to_process[[file_name]]$cell_line == "CLB_Ma") {
    clb_ma_control_samples <- intersect(sample_size_factors$sample[sample_size_factors$cell_line == "CLB_Ma" & sample_size_factors$timepoint == "control"], sample_cols)
    clb_ma_24h_samples <- intersect(sample_size_factors$sample[sample_size_factors$cell_line == "CLB_Ma" & sample_size_factors$timepoint == "24h"], sample_cols)
    
    if (length(clb_ma_control_samples) > 0 && length(clb_ma_24h_samples) > 0) {
      raw_control_mean <- rowMeans(atac_matrix_raw[, clb_ma_control_samples, drop = FALSE], na.rm = TRUE)
      raw_24h_mean <- rowMeans(atac_matrix_raw[, clb_ma_24h_samples, drop = FALSE], na.rm = TRUE)
      raw_fc <- log2((raw_24h_mean + 0.001) / (raw_control_mean + 0.001))
      
      norm_control_mean <- rowMeans(atac_matrix_normalized[, clb_ma_control_samples, drop = FALSE], na.rm = TRUE)
      norm_24h_mean <- rowMeans(atac_matrix_normalized[, clb_ma_24h_samples, drop = FALSE], na.rm = TRUE)
      norm_fc <- log2((norm_24h_mean + 0.001) / (norm_control_mean + 0.001))
      
      cat(sprintf("  Raw data - Median log2FC: %.3f\n", median(raw_fc, na.rm = TRUE)))
      cat(sprintf("  Normalized - Median log2FC: %.3f\n", median(norm_fc, na.rm = TRUE)))
    }
  }
  
  normalized_list_of_SEs[[file_name]] <- atac_scores_normalized
  cat(sprintf("  Completed normalization for %d regions\n", nrow(atac_scores_normalized)))
}

cat("\nGenome-wide normalization completed for all datasets.\n")

# ============================================================================
# STATISTICAL ANALYSIS AND FOLD CHANGE CALCULATIONS
# ============================================================================

# Function to calculate log2 fold changes
calc_log2fc <- function(matrix, control_cols, treatment_cols) {
  control_mean <- rowMeans(matrix[, control_cols, drop = FALSE], na.rm = TRUE)
  treatment_mean <- rowMeans(matrix[, treatment_cols, drop = FALSE], na.rm = TRUE)
  log2((treatment_mean + 0.001) / (control_mean + 0.001))
}

# Function to process a single dataset
process_dataset <- function(data, file_name, file_info) {
  cat(sprintf("\nProcessing %s...\n", file_info$description))
  
  # Get sample columns (excluding coordinate and metadata columns)
  sample_cols <- colnames(data)[4:(ncol(data)-3)]
  
  # Extract sample information from column names
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
  
  # Create accessibility matrix (regions x samples)
  accessibility_matrix <- as.matrix(data[, sample_cols])
  rownames(accessibility_matrix) <- paste0(file_name, "_SE_", 1:nrow(accessibility_matrix))
  
  # Calculate means for each condition
  condition_means <- sample_info %>%
    group_by(cell_line, timepoint) %>%
    summarise(
      samples = list(sample),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      mean_accessibility = list(rowMeans(accessibility_matrix[, samples], na.rm = TRUE))
    ) %>%
    ungroup()
  
  # Create condition matrix
  condition_matrix <- do.call(cbind, lapply(1:nrow(condition_means), function(i) {
    condition_means$mean_accessibility[[i]]
  }))
  colnames(condition_matrix) <- paste(condition_means$cell_line, condition_means$timepoint, sep = "_")
  rownames(condition_matrix) <- rownames(accessibility_matrix)
  
  # Calculate fold changes for each cell line
  results_list <- list()
  
  for (cell_line in unique(sample_info$cell_line)) {
    control_col <- paste0(cell_line, "_control")
    h24_col <- paste0(cell_line, "_24h") 
    h48_col <- paste0(cell_line, "_48h")
    
    if (all(c(control_col, h24_col, h48_col) %in% colnames(condition_matrix))) {
      
      fc_24h <- calc_log2fc(condition_matrix, control_col, h24_col)
      fc_48h <- calc_log2fc(condition_matrix, control_col, h48_col)
      
      results_list[[cell_line]] <- data.frame(
        se_id = rownames(condition_matrix),
        chr = data$chr,
        start = data$start,
        end = data$end,
        control_accessibility = condition_matrix[, control_col],
        h24_accessibility = condition_matrix[, h24_col],
        h48_accessibility = condition_matrix[, h48_col],
        log2fc_24h = fc_24h,
        log2fc_48h = fc_48h,
        cell_line = cell_line,
        file_type = file_name,
        phenotype = file_info$phenotype,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine results for this dataset
  dataset_results <- do.call(rbind, results_list)
  
  # Define significant changes (>1.5-fold change)
  log2_threshold <- log2(1.5)
  
  # Classify changes
  dataset_results$change_24h <- case_when(
    dataset_results$log2fc_24h > log2_threshold ~ "Increased",
    dataset_results$log2fc_24h < -log2_threshold ~ "Decreased", 
    TRUE ~ "No change"
  )
  
  dataset_results$change_48h <- case_when(
    dataset_results$log2fc_48h > log2_threshold ~ "Increased",
    dataset_results$log2fc_48h < -log2_threshold ~ "Decreased",
    TRUE ~ "No change"
  )
  
  # Handle NaN values
  dataset_results <- dataset_results %>%
    mutate(across(where(is.numeric), ~ifelse(is.nan(.), 0, .)))
  
  return(dataset_results)
}

# Process all datasets
all_dataset_results <- list()

for (file_name in names(normalized_list_of_SEs)) {
  file_info <- files_to_process[[file_name]]
  dataset_result <- process_dataset(normalized_list_of_SEs[[file_name]], file_name, file_info)
  all_dataset_results[[file_name]] <- dataset_result
}

# Combine all results
all_results <- do.call(rbind, all_dataset_results)
rownames(all_results) <- NULL

# Statistical summary
cat("\n=== UNIQUE SUPER ENHANCER ACCESSIBILITY CHANGES ===\n")

for (file_name in names(all_dataset_results)) {
  file_info <- files_to_process[[file_name]]
  cat(sprintf("\n%s:\n", file_info$description))
  
  dataset_data <- all_dataset_results[[file_name]]
  
  for (cell_line in unique(dataset_data$cell_line)) {
    cell_data <- dataset_data[dataset_data$cell_line == cell_line, ]
    if (nrow(cell_data) == 0) next
    
    cat(sprintf("\n  %s:\n", cell_line))
    cat("    24h treatment:\n")
    print(table(cell_data$change_24h))
    cat("    48h treatment:\n") 
    print(table(cell_data$change_48h))
    
    # Median fold changes
    cat(sprintf("    Median log2FC 24h: %.3f\n", median(cell_data$log2fc_24h, na.rm = TRUE)))
    cat(sprintf("    Median log2FC 48h: %.3f\n", median(cell_data$log2fc_48h, na.rm = TRUE)))
  }
}

cat(sprintf("\nTotal unique super enhancers analyzed: %d\n", nrow(all_results)))
cat(sprintf("Datasets: %s\n", paste(names(all_dataset_results), collapse = ", ")))

# Create summary statistics
summary_stats <- all_results %>%
  group_by(file_type, phenotype, cell_line) %>%
  summarise(
    total_ses = n(),
    decreased_24h = sum(change_24h == "Decreased"),
    increased_24h = sum(change_24h == "Increased"),
    decreased_48h = sum(change_48h == "Decreased"),
    increased_48h = sum(change_48h == "Increased"),
    median_fc_24h = median(log2fc_24h, na.rm = TRUE),
    median_fc_48h = median(log2fc_48h, na.rm = TRUE),
    mean_control_accessibility = mean(control_accessibility, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

library(gridExtra)
library(textshape)

# Function to create comparative boxplots for each dataset type
create_dataset_boxplots <- function() {
  cat("\nCreating comparative boxplots...\n")
  
  # Prepare data for visualization
  plot_data <- all_results %>%
    select(se_id, cell_line, phenotype, file_type, control_accessibility, h24_accessibility, h48_accessibility) %>%
    pivot_longer(cols = c(control_accessibility, h24_accessibility, h48_accessibility),
                 names_to = "timepoint", values_to = "accessibility") %>%
    mutate(
      timepoint_clean = case_when(
        timepoint == "control_accessibility" ~ "Control",
        timepoint == "h24_accessibility" ~ "24h",
        timepoint == "h48_accessibility" ~ "48h"
      ),
      timepoint_clean = factor(timepoint_clean, levels = c("Control", "24h", "48h")),
      cell_phenotype = paste(cell_line, phenotype, sep = "_"),
      accessibility = ifelse(is.nan(accessibility), 0, accessibility)
    )
  
  # Create boxplot comparing all unique SE types
  p_all <- ggplot(plot_data, aes(x = timepoint_clean, y = accessibility, fill = cell_phenotype)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.3) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("Control", "24h"), c("Control", "48h"), c("24h", "48h")),
      label = "p.signif",
      size = 3
    ) +
    facet_wrap(~cell_phenotype, scales = "free_y", ncol = 2) +
    scale_fill_viridis_d(option = "plasma", name = "SE Type") +
    labs(
      x = "Treatment Duration",
      y = "ATAC-seq Accessibility Score",
      title = "Unique Super Enhancer Accessibility Changes",
      subtitle = "K975 Treatment Effects on Cell-Type Specific Super Enhancers"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  print(p_all)
  ggsave(file.path(results_dir, "unique_SEs_all_boxplot.pdf"), p_all, width = 12, height = 10)
  
  return(p_all)
}

# Function to create heatmaps for each dataset type
create_dataset_heatmaps <- function() {
  cat("\nCreating heatmaps...\n")
  
  for (file_name in names(all_dataset_results)) {
    file_info <- files_to_process[[file_name]]
    dataset_results <- all_dataset_results[[file_name]]
    
    # Prepare matrix for heatmap (log2 fold changes)
    heatmap_matrix <- dataset_results %>%
      select(se_id, cell_line, log2fc_24h, log2fc_48h) %>%
      pivot_longer(cols = c(log2fc_24h, log2fc_48h), names_to = "timepoint", values_to = "log2fc") %>%
      unite("condition", cell_line, timepoint, sep = "_") %>%
      pivot_wider(names_from = condition, values_from = log2fc) %>%
      column_to_rownames("se_id") %>%
      as.matrix()
    
    # Skip if no data
    if (nrow(heatmap_matrix) == 0) next
    
    # Create heatmap
    pdf(file.path(results_dir, paste0("heatmap_", file_name, ".pdf")), width = 8, height = max(6, nrow(heatmap_matrix)/10))
    pheatmap(
      heatmap_matrix,
      scale = "none",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      breaks = seq(-2, 2, length.out = 101),
      main = paste("Accessibility Changes:", file_info$description, "\n(Log2 Fold Change vs Control)"),
      fontsize = 10,
      angle_col = 45
    )
    dev.off()
    
    cat(sprintf("  Created heatmap for %s\n", file_info$description))
  }
}

# Function to create scatter plots comparing 24h vs 48h changes
create_scatter_plots <- function() {
  cat("\nCreating scatter plots...\n")
  
  # Scatter plot comparing 24h vs 48h changes for all datasets
  p_scatter <- ggplot(all_results, aes(x = log2fc_24h, y = log2fc_48h, color = paste(cell_line, phenotype, sep = "_"))) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    facet_wrap(~paste(cell_line, phenotype, sep = "_"), scales = "free") +
    scale_color_viridis_d(option = "plasma", name = "SE Type") +
    labs(
      x = "Log2FC Accessibility (24h vs Control)",
      y = "Log2FC Accessibility (48h vs Control)", 
      title = "Unique Super Enhancer Accessibility Changes Over Time",
      subtitle = "K975 Treatment Effect Comparison"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  print(p_scatter)
  ggsave(file.path(results_dir, "unique_SEs_scatter_plot.pdf"), p_scatter, width = 12, height = 8)
  
  return(p_scatter)
}

# Function to create log2FC distribution boxplots
create_log2fc_boxplots <- function() {
  cat("\nCreating log2FC distribution plots...\n")
  
  # Prepare data for log2FC boxplots
  log2fc_data <- all_results %>%
    select(se_id, cell_line, phenotype, file_type, log2fc_24h, log2fc_48h) %>%
    pivot_longer(cols = c(log2fc_24h, log2fc_48h), names_to = "timepoint", values_to = "log2fc") %>%
    mutate(
      timepoint = gsub("log2fc_", "", timepoint),
      cell_phenotype = paste(cell_line, phenotype, sep = "_")
    )
  
  log2_threshold <- log2(1.5)
  
  p_log2fc <- ggplot(log2fc_data, aes(x = timepoint, y = log2fc, fill = timepoint)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dotted", alpha = 0.5) +
    facet_wrap(~cell_phenotype, scales = "free_y") +
    scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.8) +
    labs(
      x = "Treatment Duration",
      y = "Log2 Fold Change in Accessibility",
      title = "Distribution of Unique Super Enhancer Accessibility Changes",
      subtitle = "Dashed lines: ±1.5-fold change threshold"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "none"
    )
  
  print(p_log2fc)
  ggsave(file.path(results_dir, "unique_SEs_log2fc_boxplot.pdf"), p_log2fc, width = 12, height = 8)
  
  return(p_log2fc)
}

# Function to create summary comparison plot
create_summary_comparison <- function() {
  cat("\nCreating summary comparison plot...\n")
  
  # Prepare summary data for plotting
  summary_plot_data <- summary_stats %>%
    select(file_type, phenotype, cell_line, increased_24h, decreased_24h, increased_48h, decreased_48h) %>%
    pivot_longer(cols = c(increased_24h, decreased_24h, increased_48h, decreased_48h),
                 names_to = "change_type", values_to = "count") %>%
    separate(change_type, into = c("direction", "timepoint"), sep = "_") %>%
    mutate(
      direction = factor(direction, levels = c("increased", "decreased")),
      timepoint = paste0(timepoint, "h"),
      cell_phenotype = paste(cell_line, phenotype, sep = "_")
    )
  
  p_summary <- ggplot(summary_plot_data, aes(x = timepoint, y = count, fill = direction)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~cell_phenotype, scales = "free_y") +
    scale_fill_manual(values = c("increased" = "#E31A1C", "decreased" = "#1F78B4"), name = "Change Direction") +
    labs(
      x = "Treatment Duration",
      y = "Number of Super Enhancers",
      title = "Summary of Unique Super Enhancer Changes",
      subtitle = "Count of significantly changed regions (>1.5-fold)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p_summary)
  ggsave(file.path(results_dir, "unique_SEs_summary_comparison.pdf"), p_summary, width = 12, height = 8)
  
  return(p_summary)
}

# Execute all visualizations
cat("\n=== CREATING VISUALIZATIONS ===\n")
plot_boxplot <- create_dataset_boxplots()
create_dataset_heatmaps()
plot_scatter <- create_scatter_plots()
plot_log2fc <- create_log2fc_boxplots()
plot_summary <- create_summary_comparison()

cat("\nAll visualizations completed and saved to:", results_dir, "\n")

# ============================================================================
# RESULTS EXPORT
# ============================================================================

cat("\n=== EXPORTING RESULTS ===\n")

# Export detailed results for all datasets
write.csv(all_results, file.path(results_dir, "unique_SEs_detailed_results.csv"), row.names = FALSE)
cat("Exported detailed results for all unique super enhancers\n")

# Export summary statistics
write.csv(summary_stats, file.path(results_dir, "unique_SEs_summary_statistics.csv"), row.names = FALSE)
cat("Exported summary statistics\n")

# Export individual dataset results
for (file_name in names(all_dataset_results)) {
  file_info <- files_to_process[[file_name]]
  dataset_results <- all_dataset_results[[file_name]]
  
  output_file <- file.path(results_dir, paste0("detailed_results_", file_name, ".csv"))
  write.csv(dataset_results, output_file, row.names = FALSE)
  cat(sprintf("Exported %s results (%d regions)\n", file_info$description, nrow(dataset_results)))
}

# Create and export top affected regions for each dataset
top_affected_all <- list()

for (file_name in names(all_dataset_results)) {
  file_info <- files_to_process[[file_name]]
  dataset_results <- all_dataset_results[[file_name]]
  
  # Get top 50 most affected regions
  top_affected <- dataset_results %>%
    mutate(
      max_abs_change = pmax(abs(log2fc_24h), abs(log2fc_48h), na.rm = TRUE),
      primary_effect = ifelse(abs(log2fc_24h) > abs(log2fc_48h), "24h", "48h")
    ) %>%
    arrange(desc(max_abs_change)) %>%
    head(50) %>%
    mutate(dataset = file_name, phenotype = file_info$phenotype)
  
  top_affected_all[[file_name]] <- top_affected
  
  # Export individual top affected
  output_file <- file.path(results_dir, paste0("top50_affected_", file_name, ".csv"))
  write.csv(top_affected, output_file, row.names = FALSE)
  cat(sprintf("Exported top 50 affected regions for %s\n", file_info$description))
}

# Combine and export all top affected regions
all_top_affected <- do.call(rbind, top_affected_all)
write.csv(all_top_affected, file.path(results_dir, "top50_affected_all_datasets.csv"), row.names = FALSE)
cat("Exported combined top 50 affected regions from all datasets\n")

# Create analysis metadata file
metadata <- list(
  analysis_date = Sys.Date(),
  script_name = "k975_ATACseq_compare_unique_SEs.R",
  datasets_analyzed = names(files_to_process),
  total_regions = nrow(all_results),
  cell_lines = unique(all_results$cell_line),
  phenotypes = unique(all_results$phenotype),
  treatment_timepoints = c("control", "24h", "48h"),
  log2fc_threshold = log2(1.5),
  normalization_method = "genome-wide size factors",
  output_directory = results_dir
)

# Write metadata as JSON-like format
metadata_file <- file.path(results_dir, "analysis_metadata.txt")
writeLines(c(
  "=== K975 ATAC-seq Unique Super Enhancer Analysis Metadata ===",
  "",
  paste("Analysis Date:", metadata$analysis_date),
  paste("Script:", metadata$script_name),
  paste("Total Regions Analyzed:", metadata$total_regions),
  paste("Cell Lines:", paste(metadata$cell_lines, collapse = ", ")),
  paste("Phenotypes:", paste(metadata$phenotypes, collapse = ", ")),
  paste("Treatment Timepoints:", paste(metadata$treatment_timepoints, collapse = ", ")),
  paste("Log2FC Threshold:", metadata$log2fc_threshold),
  paste("Normalization Method:", metadata$normalization_method),
  paste("Output Directory:", metadata$output_directory),
  "",
  "Datasets Analyzed:",
  sapply(names(files_to_process), function(x) {
    paste("  -", x, ":", files_to_process[[x]]$description)
  }),
  "",
  "Output Files Generated:",
  "  - unique_SEs_detailed_results.csv",
  "  - unique_SEs_summary_statistics.csv", 
  "  - detailed_results_[dataset].csv (for each dataset)",
  "  - top50_affected_[dataset].csv (for each dataset)",
  "  - top50_affected_all_datasets.csv",
  "  - unique_SEs_all_boxplot.pdf",
  "  - heatmap_[dataset].pdf (for each dataset)",
  "  - unique_SEs_scatter_plot.pdf",
  "  - unique_SEs_log2fc_boxplot.pdf",
  "  - unique_SEs_summary_comparison.pdf"
), metadata_file)

cat("Exported analysis metadata\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", results_dir))
cat(sprintf("Total unique super enhancers analyzed: %d\n", nrow(all_results)))
cat(sprintf("Datasets: %s\n", paste(names(files_to_process), collapse = ", ")))

# Print final summary table
cat("\n=== FINAL SUMMARY BY DATASET ===\n")
final_summary <- all_results %>%
  group_by(file_type, phenotype, cell_line) %>%
  summarise(
    n_regions = n(),
    regions_increased_24h = sum(change_24h == "Increased"),
    regions_decreased_24h = sum(change_24h == "Decreased"),
    regions_increased_48h = sum(change_48h == "Increased"),
    regions_decreased_48h = sum(change_48h == "Decreased"),
    median_log2fc_24h = round(median(log2fc_24h, na.rm = TRUE), 3),
    median_log2fc_48h = round(median(log2fc_48h, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(file_type, cell_line)

print(final_summary)

# Top changing regions summary
cat("\n=== TOP 5 MOST AFFECTED REGIONS OVERALL ===\n")
top5_overall <- all_results %>%
  mutate(max_abs_change = pmax(abs(log2fc_24h), abs(log2fc_48h), na.rm = TRUE)) %>%
  arrange(desc(max_abs_change)) %>%
  head(5) %>%
  select(se_id, chr, start, end, cell_line, phenotype, log2fc_24h, log2fc_48h, max_abs_change)

print(top5_overall)

cat("\nAnalysis completed successfully!\n")
cat("All results, visualizations, and metadata have been exported.\n")