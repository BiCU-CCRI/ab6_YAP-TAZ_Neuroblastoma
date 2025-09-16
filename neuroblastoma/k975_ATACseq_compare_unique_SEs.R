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


# Load raw ATAC-seq accessibility scores at super enhancers
files_to_process <- list(
  SKN_A_u = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_skn_ADRu.tab",
  SKN_M_u = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_skn_MESu.tab",
  
  CLMB_A_u = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_clbm_ADRu.tab",
  CLMB_M_u = "~/workspace/neuroblastoma/data/ATACseq/bigWigs/scores_per_transcript_clbm_MESu.tab"
)
cat("Loading raw ATAC-seq accessibility scores...\n")
list_of_SEs <- list()

for (file_path in names(files_to_process)){
   atac_scores_raw <- read.table(
    file.path(files_to_process[[file_path]]),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
   colnames(atac_scores_raw) <- gsub("'", "", colnames(atac_scores_raw))
   colnames(atac_scores_raw) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores_raw))
   
  list_of_SEs[[file_path]] <- atac_scores_raw
}

# Clean column names and extract metadata
cat(sprintf("Loaded raw accessibility data for %d super enhancer regions\n", nrow(atac_scores_raw)))








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

# Apply normalization  super enhancer data

sample_cols <- colnames(atac_scores_raw)[4:ncol(atac_scores_raw)]
atac_matrix_raw <- as.matrix(atac_scores_raw[, sample_cols])

# Create normalized matrix
atac_matrix_normalized <- atac_matrix_raw
for (i in seq_along(sample_cols)) {
  sample_name <- sample_cols[i]
  size_factor <- sample_size_factors$size_factor[sample_size_factors$sample == sample_name]
  
  if (length(size_factor) == 1 && !is.na(size_factor)) {
    atac_matrix_normalized[, i] <- atac_matrix_raw[, i] / size_factor
    cat(sprintf("Applied size factor %.3f to %s\n", size_factor, sample_name))
  }
}

# Create final normalized data frame
atac_scores <- data.frame(
  atac_scores_raw[, 1:3],  # Keep coordinate columns
  atac_matrix_normalized
)

# Create genomic coordinates for super enhancers
atac_coords <- GRanges(
  seqnames = atac_scores$chr,
  ranges = IRanges(start = atac_scores$start, end = atac_scores$end),
  se_id = paste0("SE_", 1:nrow(atac_scores))
)

# Validate normalization effect
clb_ma_control_samples <- sample_size_factors$sample[sample_size_factors$cell_line == "CLB_Ma" & sample_size_factors$timepoint == "control"]
clb_ma_24h_samples <- sample_size_factors$sample[sample_size_factors$cell_line == "CLB_Ma" & sample_size_factors$timepoint == "24h"]

if (length(clb_ma_control_samples) > 0 && length(clb_ma_24h_samples) > 0) {
  # Calculate fold changes before and after normalization
  raw_control_mean <- rowMeans(atac_matrix_raw[, clb_ma_control_samples], na.rm = TRUE)
  raw_24h_mean <- rowMeans(atac_matrix_raw[, clb_ma_24h_samples], na.rm = TRUE)
  raw_fc <- log2((raw_24h_mean + 0.001) / (raw_control_mean + 0.001))
  
  norm_control_mean <- rowMeans(atac_matrix_normalized[, clb_ma_control_samples], na.rm = TRUE)
  norm_24h_mean <- rowMeans(atac_matrix_normalized[, clb_ma_24h_samples], na.rm = TRUE)
  norm_fc <- log2((norm_24h_mean + 0.001) / (norm_control_mean + 0.001))
  
  cat("\nNormalization validation for CLB-Ma 24h vs control:\n")
  cat(sprintf("Raw data - Median log2FC: %.3f, %% regions >0.5 log2FC: %.1f%%\n",
              median(raw_fc, na.rm = TRUE),
              100 * sum(raw_fc > 0.5, na.rm = TRUE) / length(raw_fc)))
  cat(sprintf("Normalized - Median log2FC: %.3f, %% regions >0.5 log2FC: %.1f%%\n",
              median(norm_fc, na.rm = TRUE),
              100 * sum(norm_fc > 0.5, na.rm = TRUE) / length(norm_fc)))
}

cat(sprintf("\nFinal normalized accessibility data for %d ADR super enhancer regions\n", nrow(atac_scores)))