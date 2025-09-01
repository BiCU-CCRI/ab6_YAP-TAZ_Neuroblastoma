# K975 ATAC-seq Super Enhancer Accessibility Analysis
# Analysis of accessibility changes at MES super enhancers after K975 treatment
# Comparing control vs 24h vs 48h treatment in CM (CLB-Ma) and SH (SK-N-SH) cell lines

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
results_dir <- "~/workspace/neuroblastoma/results/ATAC-seq_drug_treatment"

# Create results directory if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# 1. Load ATAC-seq accessibility scores at super enhancers
cat("Loading ATAC-seq accessibility scores...\n")
atac_scores <- read.table(
  file.path(data_dir, "ATACseq/bigWigs/scores_per_transcript.tab"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Clean column names and extract metadata
colnames(atac_scores) <- gsub("'", "", colnames(atac_scores))
colnames(atac_scores) <- gsub("\\.mLb\\.clN\\.bigWig", "", colnames(atac_scores))

# Create genomic coordinates for super enhancers
atac_coords <- GRanges(
  seqnames = atac_scores$chr,
  ranges = IRanges(start = atac_scores$start, end = atac_scores$end),
  se_id = paste0("SE_", 1:nrow(atac_scores))
)

cat(sprintf("Loaded accessibility data for %d super enhancer regions\n", nrow(atac_scores)))

# 2. Load MES super enhancer annotations from HOMER annotation file
cat("Loading MES super enhancer annotations...\n")
mes_se_file <- "~/workspace/neuroblastoma/temp_results/BEDs/CLB_SKN_M_homer_annot.csv"

if (file.exists(mes_se_file)) {
  # Read HOMER annotation file
  mes_se <- read.csv(mes_se_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract coordinates (remove 'chr' prefix if present for consistency)
  mes_se$Chr <- gsub("chr", "", mes_se$Chr)
  
  # Create GRanges for MES super enhancers
  mes_se_gr <- GRanges(
    seqnames = mes_se$Chr,
    ranges = IRanges(start = mes_se$Start, end = mes_se$End),
    peak_id = mes_se[,1],  # First column contains PeakID
    annotation = mes_se$Annotation,
    gene_name = mes_se$Gene.Name,
    distance_to_tss = mes_se$Distance.to.TSS
  )
  
  cat(sprintf("Loaded %d MES super enhancers with HOMER annotations\n", length(mes_se_gr)))
  cat(sprintf("Annotation types: %s\n", paste(unique(mes_se$Annotation)[1:min(5, length(unique(mes_se$Annotation)))], collapse = ", ")))
  
  # COORDINATE VALIDATION CHECK
  cat("\n=== COORDINATE VALIDATION ===\n")
  
  # Create comparable coordinate strings for both datasets
  atac_coord_strings <- paste(atac_scores$chr, atac_scores$start, atac_scores$end, sep = ":")
  homer_coord_strings <- paste(mes_se$Chr, mes_se$Start, mes_se$End, sep = ":")
  
  # Check for exact matches
  exact_matches <- intersect(atac_coord_strings, homer_coord_strings)
  cat(sprintf("Exact coordinate matches: %d out of %d ATAC regions (%d HOMER regions)\n", 
              length(exact_matches), length(atac_coord_strings), length(homer_coord_strings)))
  
  # Check for overlaps using GenomicRanges
  overlaps <- findOverlaps(atac_coords, mes_se_gr)
  overlap_count <- length(unique(queryHits(overlaps)))
  cat(sprintf("Overlapping regions: %d out of %d ATAC regions\n", overlap_count, length(atac_coords)))
  
  # Detailed validation
  if (length(exact_matches) == 0 && overlap_count == 0) {
    cat("❌ ERROR: No coordinate matches found between ATAC data and HOMER annotations!\n")
    cat("This suggests the coordinate systems may be different.\n")
    
    # Diagnostic information
    cat("\nDiagnostic information:\n")
    cat("ATAC coordinates (first 5):\n")
    print(head(atac_coord_strings, 5))
    cat("\nHOMER coordinates (first 5):\n")
    print(head(homer_coord_strings, 5))
    
    # Check if chromosome naming is the issue
    atac_chrs <- unique(atac_scores$chr)
    homer_chrs <- unique(mes_se$Chr)
    cat(sprintf("\nATAC chromosomes: %s\n", paste(head(atac_chrs, 10), collapse = ", ")))
    cat(sprintf("HOMER chromosomes: %s\n", paste(head(homer_chrs, 10), collapse = ", ")))
    
    stop("Coordinate validation failed. Please check that ATAC-seq and HOMER annotation files refer to the same genomic regions.")
    
  } else if (length(exact_matches) > 0) {
    cat("✅ VALIDATION PASSED: Found exact coordinate matches\n")
    cat(sprintf("Match rate: %.1f%% of ATAC regions have exact matches in HOMER annotations\n", 
                100 * length(exact_matches) / length(atac_coord_strings)))
    
    # Create mapping for matched regions
    atac_matched_idx <- match(exact_matches, atac_coord_strings)
    homer_matched_idx <- match(exact_matches, homer_coord_strings)
    
    coordinate_mapping <- data.frame(
      atac_index = atac_matched_idx,
      homer_index = homer_matched_idx,
      coordinates = exact_matches,
      stringsAsFactors = FALSE
    )
    
  } else if (overlap_count > 0) {
    cat("⚠️  PARTIAL MATCH: Found overlapping regions but no exact matches\n")
    cat("This may indicate slightly different coordinate systems or region boundaries.\n")
    cat(sprintf("Overlap rate: %.1f%% of ATAC regions overlap with HOMER annotations\n",
                100 * overlap_count / length(atac_coords)))
    
    # Create mapping for overlapping regions
    overlaps_df <- as.data.frame(overlaps)
    coordinate_mapping <- data.frame(
      atac_index = overlaps_df$queryHits,
      homer_index = overlaps_df$subjectHits,
      coordinates = paste(atac_coord_strings[overlaps_df$queryHits], 
                         homer_coord_strings[overlaps_df$subjectHits], sep = " <-> "),
      stringsAsFactors = FALSE
    )
  }
  
  # Save coordinate mapping for reference
  if (exists("coordinate_mapping")) {
    write.csv(coordinate_mapping, file.path(results_dir, "coordinate_mapping_validation.csv"), row.names = FALSE)
    cat("Coordinate mapping saved to coordinate_mapping_validation.csv\n")
  }
  
} else {
  cat("MES super enhancer annotation file not found, proceeding with all regions as potential SEs\n")
  mes_se_gr <- atac_coords  # Use all regions if SE annotation not found
  coordinate_mapping <- NULL
}

# 3. Prepare data matrix and metadata
cat("Preparing data matrix...\n")

# Extract sample information from column names
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

# Create accessibility matrix (regions x samples)
accessibility_matrix <- as.matrix(atac_scores[, sample_cols])
rownames(accessibility_matrix) <- paste0("SE_", 1:nrow(accessibility_matrix))

# 4. Calculate mean accessibility per condition
cat("Calculating condition means...\n")

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

# 5. Calculate fold changes relative to control
cat("Calculating fold changes...\n")

# Function to calculate log2 fold changes
calc_log2fc <- function(matrix, control_cols, treatment_cols) {
  control_mean <- rowMeans(matrix[, control_cols, drop = FALSE], na.rm = TRUE)
  treatment_mean <- rowMeans(matrix[, treatment_cols, drop = FALSE], na.rm = TRUE)
  log2((treatment_mean + 0.001) / (control_mean + 0.001))
}

# Calculate fold changes for each cell line
results_list <- list()

for (cell_line in c("CLB_Ma", "SK_N_SH")) {
  control_col <- paste0(cell_line, "_control")
  h24_col <- paste0(cell_line, "_24h") 
  h48_col <- paste0(cell_line, "_48h")
  
  if (all(c(control_col, h24_col, h48_col) %in% colnames(condition_matrix))) {
    
    fc_24h <- calc_log2fc(condition_matrix, control_col, h24_col)
    fc_48h <- calc_log2fc(condition_matrix, control_col, h48_col)
    
    results_list[[cell_line]] <- data.frame(
      se_id = rownames(condition_matrix),
      chr = atac_scores$chr,
      start = atac_scores$start,
      end = atac_scores$end,
      control_accessibility = condition_matrix[, control_col],
      h24_accessibility = condition_matrix[, h24_col],
      h48_accessibility = condition_matrix[, h48_col],
      log2fc_24h = fc_24h,
      log2fc_48h = fc_48h,
      cell_line = cell_line,
      stringsAsFactors = FALSE
    )
  }
}

# Combine results
all_results <- do.call(rbind, results_list)

# 6. Statistical analysis and visualization
cat("Performing statistical analysis...\n")

# Define significant changes (>1.5-fold change)
log2_threshold <- log2(1.5)

# Classify changes
all_results$change_24h <- case_when(
  all_results$log2fc_24h > log2_threshold ~ "Increased",
  all_results$log2fc_24h < -log2_threshold ~ "Decreased", 
  TRUE ~ "No change"
)

all_results$change_48h <- case_when(
  all_results$log2fc_48h > log2_threshold ~ "Increased",
  all_results$log2fc_48h < -log2_threshold ~ "Decreased",
  TRUE ~ "No change"
)

# Summary statistics
cat("=== SUPER ENHANCER ACCESSIBILITY CHANGES ===\n")
for (cell_line in c("CLB_Ma", "SK_N_SH")) {
  cat(sprintf("\n%s (%s):\n", cell_line, ifelse(cell_line == "CLB_Ma", "CM", "SH")))
  
  cell_data <- all_results[all_results$cell_line == cell_line, ]
  
  cat("24h treatment:\n")
  print(table(cell_data$change_24h))
  
  cat("48h treatment:\n") 
  print(table(cell_data$change_48h))
  
  # Median fold changes
  cat(sprintf("Median log2FC 24h: %.3f\n", median(cell_data$log2fc_24h, na.rm = TRUE)))
  cat(sprintf("Median log2FC 48h: %.3f\n", median(cell_data$log2fc_48h, na.rm = TRUE)))
}

# 7. Create visualizations

# Heatmap of accessibility changes
cat("Creating visualizations...\n")

# Prepare matrix for heatmap (log2 fold changes)
heatmap_matrix <- all_results %>%
  dplyr::select(se_id, cell_line, log2fc_24h, log2fc_48h) %>%
  pivot_longer(cols = c(log2fc_24h, log2fc_48h), names_to = "timepoint", values_to = "log2fc") %>%
  unite("condition", cell_line, timepoint, sep = "_") %>%
  pivot_wider(names_from = condition, values_from = log2fc) %>%
  textshape::column_to_rownames("se_id") %>%
  as.matrix()

# Create heatmap
pdf(file.path(results_dir, "K975_SE_accessibility_heatmap.pdf"), width = 8, height = 12)
pheatmap(
  heatmap_matrix,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  #color = colorRampPalette(rev(RdBu(100)))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "Super Enhancer Accessibility Changes\n(Log2 Fold Change vs Control)",
  fontsize = 10,
  angle_col = 45
)
dev.off()

# Scatter plot comparing 24h vs 48h changes
p1 <- ggplot(all_results, aes(x = log2fc_24h, y = log2fc_48h)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~cell_line, labeller = as_labeller(c("CLB_Ma" = "CLB-Ma (CM)", "SK_N_SH" = "SK-N-SH (SH)"))) +
  labs(
    x = "Log2FC Accessibility (24h vs Control)",
    y = "Log2FC Accessibility (48h vs Control)", 
    title = "Super Enhancer Accessibility Changes Over Time",
    subtitle = "K975 Treatment Effect"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(file.path(results_dir, "K975_SE_accessibility_scatter.pdf"), p1, width = 10, height = 6)

# Box plots showing distribution of changes
plot_data <- all_results %>%
  dplyr::select(se_id, cell_line, log2fc_24h, log2fc_48h) %>%
  pivot_longer(cols = c(log2fc_24h, log2fc_48h), names_to = "timepoint", values_to = "log2fc") %>%
  mutate(
    timepoint = gsub("log2fc_", "", timepoint),
    cell_line_label = ifelse(cell_line == "CLB_Ma", "CLB-Ma (CM)", "SK-N-SH (SH)")
  )

p2 <- ggplot(plot_data, aes(x = timepoint, y = log2fc, fill = timepoint)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = c(-log2_threshold, log2_threshold), linetype = "dotted", alpha = 0.5) +
  facet_wrap(~cell_line_label) +
  scale_fill_viridis_d(option = "plasma", begin = 0.3, end = 0.8) +
  labs(
    x = "Treatment Duration",
    y = "Log2 Fold Change in Accessibility",
    title = "Distribution of Super Enhancer Accessibility Changes",
    subtitle = "Dashed lines: ±1.5-fold change threshold"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(results_dir, "K975_SE_accessibility_boxplot.pdf"), p2, width = 10, height = 6)

# Time course plot for top changing SEs
top_changing_ses <- all_results %>%
  mutate(max_change = pmax(abs(log2fc_24h), abs(log2fc_48h))) %>%
  arrange(desc(max_change)) %>%
  head(20)

# Prepare data for time course plot
timecourse_data <- top_changing_ses %>%
  dplyr::select(se_id, cell_line, control_accessibility, h24_accessibility, h48_accessibility) %>%
  pivot_longer(cols = c(control_accessibility, h24_accessibility, h48_accessibility),
               names_to = "timepoint", values_to = "accessibility") %>%
  mutate(
    timepoint_num = case_when(
      timepoint == "control_accessibility" ~ 0,
      timepoint == "h24_accessibility" ~ 24,
      timepoint == "h48_accessibility" ~ 48
    ),
    cell_line_label = ifelse(cell_line == "CLB_Ma", "CLB-Ma (CM)", "SK-N-SH (SH)")
  )

p3 <- ggplot(timecourse_data, aes(x = timepoint_num, y = accessibility, group = se_id)) +
  geom_line(alpha = 0.5, color = "grey60") +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "red", size = 1.5) +
  facet_wrap(~cell_line_label) +
  labs(
    x = "Treatment Duration (hours)",
    y = "ATAC-seq Accessibility Score",
    title = "Time Course of Super Enhancer Accessibility Changes",
    subtitle = "Top 20 most variable super enhancers (individual SEs in grey, trend in red)"
  ) +
  scale_x_continuous(breaks = c(0, 24, 48)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(file.path(results_dir, "K975_SE_accessibility_timecourse.pdf"), p3, width = 10, height = 6)

# 8. Export results
cat("Exporting results...\n")

# Save detailed results
write.csv(all_results, file.path(results_dir, "K975_SE_accessibility_detailed_results.csv"), row.names = FALSE)

# Save summary statistics
summary_stats <- all_results %>%
  group_by(cell_line) %>%
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

write.csv(summary_stats, file.path(results_dir, "K975_SE_accessibility_summary.csv"), row.names = FALSE)

# Create a ranked list of most affected super enhancers
top_affected <- all_results %>%
  mutate(
    max_abs_change = pmax(abs(log2fc_24h), abs(log2fc_48h)),
    primary_effect = ifelse(abs(log2fc_24h) > abs(log2fc_48h), "24h", "48h")
  ) %>%
  arrange(desc(max_abs_change)) %>%
  head(50)

# Add HOMER annotation information to top affected SEs if available
if (exists("mes_se") && nrow(mes_se) > 0 && exists("coordinate_mapping")) {
  # Create lookup table for annotations using validated coordinate mapping
  if (!is.null(coordinate_mapping) && nrow(coordinate_mapping) > 0) {
    
    annotation_lookup <- data.frame(
      chr = mes_se$Chr[coordinate_mapping$homer_index],
      start = mes_se$Start[coordinate_mapping$homer_index],
      end = mes_se$End[coordinate_mapping$homer_index],
      peak_id = mes_se[coordinate_mapping$homer_index, 1],
      annotation = mes_se$Annotation[coordinate_mapping$homer_index],
      gene_name = mes_se$Gene.Name[coordinate_mapping$homer_index],
      distance_to_tss = mes_se$Distance.to.TSS[coordinate_mapping$homer_index],
      gene_description = if("Gene.Description" %in% colnames(mes_se)) mes_se$Gene.Description[coordinate_mapping$homer_index] else NA,
      atac_index = coordinate_mapping$atac_index,
      stringsAsFactors = FALSE
    )
    
    # Add ATAC se_id for joining
    annotation_lookup$se_id <- paste0("SE_", annotation_lookup$atac_index)
    
    # Add annotation info to top affected regions using se_id
    top_affected_annotated <- top_affected %>%
      left_join(annotation_lookup %>% select(-chr, -start, -end, -atac_index), by = "se_id") %>%
      select(se_id, chr, start, end, cell_line, log2fc_24h, log2fc_48h, max_abs_change, 
             primary_effect, peak_id, annotation, gene_name, distance_to_tss, gene_description)
    
    cat(sprintf("Successfully annotated %d out of %d top affected regions\n", 
                sum(!is.na(top_affected_annotated$annotation)), nrow(top_affected_annotated)))
    
  } else {
    cat("Warning: No valid coordinate mapping found, proceeding without annotations\n")
    top_affected_annotated <- top_affected
  }
} else if (exists("mes_se") && nrow(mes_se) > 0) {
  # Fallback: try direct coordinate matching (original approach)
  cat("Using fallback coordinate matching method...\n")
  annotation_lookup <- data.frame(
    chr = mes_se$Chr,
    start = mes_se$Start,
    end = mes_se$End,
    peak_id = mes_se[,1],
    annotation = mes_se$Annotation,
    gene_name = mes_se$Gene.Name,
    distance_to_tss = mes_se$Distance.to.TSS,
    gene_description = if("Gene.Description" %in% colnames(mes_se)) mes_se$Gene.Description else NA,
    stringsAsFactors = FALSE
  )
  
  # Add annotation info to top affected regions
  top_affected_annotated <- top_affected %>%
    left_join(annotation_lookup, by = c("chr", "start", "end")) %>%
    select(se_id, chr, start, end, cell_line, log2fc_24h, log2fc_48h, max_abs_change, 
           primary_effect, peak_id, annotation, gene_name, distance_to_tss, gene_description)
} else {
  top_affected_annotated <- top_affected
}
  
  write.csv(top_affected_annotated, file.path(results_dir, "K975_top50_affected_SEs_annotated.csv"), row.names = FALSE)
  
  # Summary by annotation type
  if (nrow(top_affected_annotated) > 0) {
    annot_summary <- top_affected_annotated %>%
      filter(!is.na(annotation)) %>%
      mutate(
        annotation_simple = case_when(
          grepl("promoter", annotation, ignore.case = TRUE) ~ "Promoter",
          grepl("intron", annotation, ignore.case = TRUE) ~ "Intronic",
          grepl("intergenic", annotation, ignore.case = TRUE) ~ "Intergenic",
          grepl("exon", annotation, ignore.case = TRUE) ~ "Exonic",
          TRUE ~ "Other"
        )
      ) %>%
      group_by(annotation_simple, cell_line) %>%
      summarise(
        count = n(),
        mean_abs_change = mean(max_abs_change, na.rm = TRUE),
        .groups = "drop"
      )
    
    write.csv(annot_summary, file.path(results_dir, "K975_annotation_type_summary.csv"), row.names = FALSE)
    cat("\nAnnotation type summary for top affected SEs:\n")
    print(annot_summary)
  }
} else {
  write.csv(top_affected, file.path(results_dir, "K975_top50_affected_SEs.csv"), row.names = FALSE)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", results_dir))
cat(sprintf("Total super enhancers analyzed: %d\n", nrow(all_results) / 2))
cat(sprintf("Cell lines: %s\n", paste(unique(all_results$cell_line), collapse = ", ")))

# Print final summary
cat("\nSUMMARY OF FINDINGS:\n")
print(summary_stats)

cat("\nTop 10 most affected super enhancers:\n")
print(head(top_affected[, c("se_id", "chr", "start", "end", "cell_line", "log2fc_24h", "log2fc_48h", "max_abs_change")], 10))

cat("\nAnalysis completed successfully!\n")