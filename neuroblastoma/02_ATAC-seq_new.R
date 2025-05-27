###
# Title: ATAC-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
# 

# Set up the environment ####
deg_dir <- "~/workspace/neuroblastoma/results/ATAC-seq_drug_treatment/"

## Loading libraries ####
library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(IRanges)
library(patchwork)
import::from(stringr, str_extract, str_detect)
import::from(GenomicRanges,GRanges)
import::from(rtracklayer, liftOver, import.chain)
import::from(openxlsx, addWorksheet, writeData, saveWorkbook)
import::from(openxlsx2, read_xlsx)
import::from(.from = TFBSTools, getMatrixSet)
import::from(chromVAR, addGCBias)
import::from(.from = RColorBrewer, brewer.pal)
import::from(.from = "~/workspace/neuroblastoma/resources/utilityScripts.R",
             "generatePCA", 
#             "extract_results_DDS",
#             "meanExprsPerGroup",
             "extract_results_DDS_HIC",
             "plotVolcano",
             "get_the_poissoon_p_val")
import::from(
  .from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"),
  "filterDatasets",
  .character_only = TRUE
)

# BugFix - JASPAR2022 cannot be pulled from the web automatically, so it need to be loaded manually
download.file(url = "https://jaspar2022.genereg.net/download/database/JASPAR2022.sqlite",
              destfile = "/home/rstudio/.cache/R/JASPAR2022.sqlite")
JASPAR2022 <-  "/home/rstudio/.cache/R/JASPAR2022.sqlite"

# Analysis part ###########################################################################
# Assembling the metadata to a DESeq2 object
# We use only filtered consensus peaks. And combined annotation from genecode, homer and reg.
tmp_dir <- "~/workspace/neuroblastoma/data/ATACseq/new/"
path_to_the_ATAC_counts_data <- file.path(tmp_dir, "consensus_peaks.mLb.clN.rds")
path_annnotation_ATAC_data <- file.path(tmp_dir, "consensus_peaks.mLb.clN.annotatePeaks.txt")
ATAC_counts_data <- readRDS(path_to_the_ATAC_counts_data)
ATAC_annotation_data <- read.csv(path_annnotation_ATAC_data, sep = "\t")

colnames(ATAC_annotation_data)[1] <- "PeakId"
ATAC_annotation_data <- ATAC_annotation_data %>%
  mutate(PeakId = factor(PeakId, levels = row.names(ATAC_counts_data))) %>%
  arrange(PeakId)

#remove Scaffolds chrs
ATAC_annotation_data <- ATAC_annotation_data[ATAC_annotation_data$Chr %in% c(1:22, "X", "Y"), ]
ATAC_counts_data <- ATAC_counts_data[rownames(ATAC_counts_data) %in% ATAC_annotation_data$PeakId, ]

#add ranges
rowRanges(ATAC_counts_data) <- GenomicRanges::GRanges(seqnames = paste0("chr", ATAC_annotation_data$Chr), 
                                              ranges = IRanges::IRanges(start = ATAC_annotation_data$Start, 
                                                                        end = ATAC_annotation_data$End
                                              )
                                              
)
mcols(ATAC_counts_data) <- DataFrame(mcols(ATAC_counts_data), ATAC_annotation_data)


# WRITE METADATA 
# ATAC_metadata_df <- data.frame(matrix(nrow = dim(ATAC_counts_data)[2]))
# ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
#                                   sample_names = colnames(ATAC_counts_data)
# )
### NOT DONE YET




# Drug treatment subset
samples_to_subset <- c(
  "CM_control_ATAC_S165169_REP2",
  "CM_control_ATAC_S165169_REP1",
  "CM_24h_ATAC_S165170_REP2",
  "CM_24h_ATAC_S165170_REP1",
  "CM_48h_ATAC_S165167_REP2",
  "CM_48h_ATAC_S165167_REP1",
  "SH_control_ATAC_S165166_REP2",
  "SH_control_ATAC_S165166_REP1",
  "SH_24h_ATAC_S165164_REP2",
  "SH_24h_ATAC_S165164_REP1",
  "SH_48h_ATAC_S165165_REP2",
  "SH_48h_ATAC_S165165_REP1"
)

drg_trtm_ATAC <- ATAC_counts_data[, colData(ATAC_counts_data)$sample %in% samples_to_subset]

meta_data <- t(as.data.frame(strsplit(drg_trtm_ATAC$sample, "_")))[,c(1,2,5)]
colnames(meta_data) <- c("cell_line", "time", "replicate")
colData(drg_trtm_ATAC) <- cbind(colData(drg_trtm_ATAC), meta_data)
colData(drg_trtm_ATAC)$time <- factor(colData(drg_trtm_ATAC)$time, levels = c("control", "24h", "48h"))
colData(drg_trtm_ATAC)$cell_line <- factor(colData(drg_trtm_ATAC)$cell_line )
design(drg_trtm_ATAC) <- as.formula("~ cell_line + time")


drg_trtm_ATAC@rowRanges@elementMetadata@listData[["PeakId"]]

#
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
RNA_SEQ_data <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES")) %>%
  select(gene_symbol, Term)
mcols(drg_trtm_ATAC)$gene_category_our_RNAseq <- "ND"
mcols(drg_trtm_ATAC)$gene_category_our_RNAseq <- case_when(
  mcols(drg_trtm_ATAC)$Gene.Name %in% RNA_SEQ_data[RNA_SEQ_data$Term == "ADRN", 1] ~ "g_ADRN",
  mcols(drg_trtm_ATAC)$Gene.Name %in% RNA_SEQ_data[RNA_SEQ_data$Term == "MES", 1] ~ "g_MES",
  TRUE ~ mcols(drg_trtm_ATAC)$gene_category_our_RNAseq
)

# Add metadata - marks peaks that are around 3kb from TSS
mcols(drg_trtm_ATAC)$is_promoter_3kb <- FALSE
mcols(drg_trtm_ATAC)$is_promoter_3kb[abs(mcols(drg_trtm_ATAC)$Distance.to.TSS) <= 3000] <- TRUE


# plot full ATAC-seq PCA
drg_trtm_ATAC <- filterDatasets(drg_trtm_ATAC,
                                abs_filt = TRUE,
                                abs_filt_samples = 2
)

design(drg_trtm_ATAC) <- as.formula("~ cell_line + time")
drg_trtm_ATAC <- estimateSizeFactors(drg_trtm_ATAC)
drg_trtm_ATAC <- DESeq(drg_trtm_ATAC)
vsd <- vst(drg_trtm_ATAC, blind = F)


pca_dar <- generatePCA(transf_object = vsd, 
                       cond_interest_varPart = c("time", "cell_line"), 
                       color_variable = "time", 
                       shape_variable = "cell_line",
                       ntop_genes = 1000) +
  ggtitle("Original dataset") 


resultsNames(drg_trtm_ATAC)
coeff_name <- "time_48h_vs_control"
cond_numerator <-  "48h"
cond_denominator <-  "control"
cond_variable <- "time"
padj_cutoff = 0.05
log2FC_cutoff = 1
sample_result <- results(drg_trtm_ATAC, contrast = c(cond_variable, cond_numerator, cond_denominator))
ATAC_dds_results_48h <- extract_results_DDS(dds_object = drg_trtm_ATAC,
                                            dds_result = sample_result,
                                            coeff_name = coeff_name,
                                            cond_numerator = cond_numerator,
                                            cond_denominator = cond_denominator,
                                            cond_variable = cond_variable,
                                            padj_cutoff = padj_cutoff,
                                            log2FC_cutoff = log2FC_cutoff)
# Write the result to the xlsx file
XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_48h$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_48h$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_48h$results_all, sheet = "results_all")
openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir,  paste0("DAR_results_time_48h_vs_control.xlsx")), overwrite = T)


resultsNames(drg_trtm_ATAC)
coeff_name <- "time_24h_vs_control"
cond_numerator <-  "24h"
cond_denominator <-  "control"
cond_variable <- "time"
padj_cutoff = 0.05
log2FC_cutoff = 1
sample_result <- results(drg_trtm_ATAC, contrast = c(cond_variable, cond_numerator, cond_denominator))
ATAC_dds_results_24h <- extract_results_DDS(dds_object = drg_trtm_ATAC,
                                            dds_result = sample_result,
                                            coeff_name = coeff_name,
                                            cond_numerator = cond_numerator,
                                            cond_denominator = cond_denominator,
                                            cond_variable = cond_variable,
                                            padj_cutoff = padj_cutoff,
                                            log2FC_cutoff = log2FC_cutoff)

XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_24h$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_24h$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results_24h$results_all, sheet = "results_all")
openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir,  paste0("DAR_results_time_24h_vs_control.xlsx")), overwrite = T)



########
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["species"]] <- "9606"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- TFBSTools::getMatrixSet(JASPAR2022, opts)

# Using different subsets of peaks to check 
subsets_to_run <- c("ALL")

sbst <-  c("ALL")
# here we check if there are any peaks that have less than 5 reads across all samples.
dim(drg_trtm_ATAC)
counts_consensus_filt <- drg_trtm_ATAC[rowSums(assay(drg_trtm_ATAC)) > 5, ]
dim(counts_consensus_filt)

print(paste("procesing: ", sbst ))

# ChromVar package uses GC content to identify background peaks that are likely to be non-functional and exclude them from further analysis. 
# This is because GC-rich regions tend to have higher nucleosome occupancy and lower DNase I hypersensitivity, 
# which are indicators of closed chromatin and reduced accessibility to transcription factors.
# By identifying and removing background peaks that are likely to be non-functional, 
# ChromVar is able to focus on the regulatory regions that are most likely to be involved in transcriptional regulation. 
# This increases the sensitivity and specificity of the analysis and helps to avoid false positive results.
# Therefore, GC content is an important parameter to consider when working with the ChromVar package.
counts_consensus_filt <- chromVAR::addGCBias(counts_consensus_filt, genome = BSgenome.Hsapiens.UCSC.hg38) 

# Having corrected for bias, we can use the matchMotifs function to identify motifs under our ATACseq peaks.
# Here we supply our RangedSummarizedExperiment of counts in peaks and the genome of interest to the matchMotifs function and use the default out of matches.

motif_matches <- motifmatchr::matchMotifs(pwms = motifsToScan, 
                                          subject = counts_consensus_filt, 
                                          genome = BSgenome.Hsapiens.UCSC.hg38, 
                                          out = "matches")
motif_matches

BiocParallel::register(BiocParallel::MulticoreParam(8, progressbar = FALSE))
set.seed(42069)

# Run chromVar ----
#Following the identification of motifs in our peaks, we can perform the summarization of ATACseq signal to motifs using the computeDeviations and the computeVariability functions.

#The function computeDeviations will use a set of background peaks for normalizing the deviation scores. 
#This computation is done internally by default and not returned – to have greater control over this step, 
#a user can run the getBackgroundPeaks function themselves and pass the result to computeDeviations under the background_peaks parameter.
#Background peaks are peaks that are similar to a peak in GC content and average accessibility
#The result from getBackgroundPeaks is a matrix of indices, where each column represents the index of the peak that is a background peak.

background_peaks <- chromVAR::getBackgroundPeaks(object = counts_consensus_filt) 
access_expectation <- chromVAR::computeExpectations(object = counts_consensus_filt)
chrom_access_deviations <- chromVAR::computeDeviations(object = counts_consensus_filt, 
                                                       annotations = motif_matches,
                                                       background_peaks = background_peaks,
                                                       expectation = access_expectation)

# we can check the correlation of samples
sample_cor <- chromVAR::getSampleCorrelation(chrom_access_deviations)
annotation_row_sampleCor <- as.data.frame(SummarizedExperiment::colData(chrom_access_deviations)) %>%
  dplyr::select(sample)
pheatmap::pheatmap(as.dist(sample_cor), 
                   annotation_row = annotation_row_sampleCor, 
                   clustering_distance_rows = as.dist(1-sample_cor), 
                   clustering_distance_cols = as.dist(1-sample_cor),  annotation_names_row = TRUE, 
                    ) 

#checkign the clusterrization using tSNE
tsne_results <- chromVAR::deviationsTsne(chrom_access_deviations, threshold = 1.5, perplexity = 8, 
                                         what = "samples", shiny = FALSE)

tsne_plots <- chromVAR::plotDeviationsTsne(chrom_access_deviations, tsne_results,
                                           sample_column = "sample", shiny = FALSE)
tsne_plots

diff_acc <- chromVAR::differentialDeviations(chrom_access_deviations, groups="time", parametric = FALSE)
motif_annot <- as.data.frame(SummarizedExperiment::rowData(chrom_access_deviations)) %>%
  tibble::rownames_to_column(var = "motif") %>%
  dplyr::select(motif, name)
diff_acc_annot <- diff_acc %>%
  tibble::rownames_to_column(var = "motif") %>%
  dplyr::left_join(., motif_annot, by = "motif")
head(diff_acc_annot) 

diff_var <- chromVAR::differentialVariability(chrom_access_deviations, "time", parametric = FALSE)
diff_var_annot <- diff_var %>%
  tibble::rownames_to_column(var = "motif") %>%
  dplyr::left_join(., motif_annot, by = "motif")
head(diff_var_annot)

devZscores <- chromVAR::deviationScores(chrom_access_deviations)
devZscores_df <- tibble::rownames_to_column(as.data.frame(devZscores), var="jaspar_id")

#The function computeVariability returns a data.frame that contains the variability (standard deviation of the z scores computed above
# across all cell/samples for a set of peaks), bootstrap confidence intervals for that variability (by resampling cells/samples), 
# and a p-value for the variability being greater than the null hypothesis of 1.
chrom_access_variability <- chromVAR::computeVariability(chrom_access_deviations)
chrom_access_variability_plot <- chromVAR::plotVariability(chrom_access_variability, use_plotly = FALSE)

ggsave(chrom_access_variability_plot, 
       filename = file.path(deg_dir, paste0("ChromVAR_chrom_access_variability_corr_plot_", sbst, ".pdf")))

BiocParallel::register(BiocParallel::SerialParam())

chrom_access_variability_ord <- chrom_access_variability[order(chrom_access_variability$p_value), ]  # ordering based on p-value

#! to store
chrom_access_variability_ord_results <- chrom_access_variability_ord %>%
  tibble::rownames_to_column(var="jaspar_id") %>%
  dplyr::left_join(., devZscores_df, by = "jaspar_id")

message("Saving motif enrichment results")
ntop <- 50
topVariable <- chrom_access_variability_ord[1:ntop, ]
topVariable <- tibble::rownames_to_column(topVariable, var="jaspar_id")

devTop <- topVariable %>%
  dplyr::left_join(., devZscores_df, by = "jaspar_id")

devTop <- devTop %>% mutate(long_name = paste0(jaspar_id, "_", name))

devToPlot <- devTop %>%
  dplyr::select(long_name, 
                CM_control_ATAC_S165169_REP1,  CM_control_ATAC_S165169_REP2,  
                CM_24h_ATAC_S165170_REP1,  CM_24h_ATAC_S165170_REP2,
                CM_48h_ATAC_S165167_REP1, CM_48h_ATAC_S165167_REP2,
                SH_control_ATAC_S165166_REP1, SH_control_ATAC_S165166_REP2,
                SH_24h_ATAC_S165164_REP1, SH_24h_ATAC_S165164_REP2,
                SH_48h_ATAC_S165165_REP1, SH_48h_ATAC_S165165_REP2
                ) %>%
  tibble::column_to_rownames(var = "long_name")

annotCol_forHeatmap <- as.data.frame(colData(drg_trtm_ATAC)) %>%
  dplyr::select(time)
# annotCol_forHeatmap_colors <- list(phenotype = c(A = "#525252", M = "#fc4e2a"))
annotCol_forHeatmap_colors <- list(phenotype = c(`24h` = "#525252", `48h` = "#fc4e2a", control = "green"))



color.scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

heatmap <- pheatmap::pheatmap(as.matrix(devToPlot),
                              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                              border_color = NA,
                              #color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
                              scale = "row",
                              cluster_cols = FALSE,
                              cluster_rows = TRUE,
                              show_colnames = TRUE,
                              show_rownames = TRUE,
                              annotation_col = annotCol_forHeatmap,
                              annotation_colors = annotCol_forHeatmap_colors,  
                              fontsize_row = 9,
                              fontsize = 9,
                              main = paste0("Heatmap of signif. ATAC-seq Drug treatment control vs 24h vs 48h. ", sbst, "- used"))

ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_drug_treatment_control_vs_24h_48h_Motifs", sbst, "- used", ".png")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_drug_treatment_control_vs_24h_48h_Motifs", sbst, "- used", ".eps")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")

#separate CM and SH
devToPlot_CM <- devTop %>%
  dplyr::select(long_name, 
                CM_control_ATAC_S165169_REP1,  CM_control_ATAC_S165169_REP2,  
                CM_24h_ATAC_S165170_REP1,  CM_24h_ATAC_S165170_REP2,
                CM_48h_ATAC_S165167_REP1, CM_48h_ATAC_S165167_REP2
  ) %>%
  tibble::column_to_rownames(var = "long_name")

heatmap <- pheatmap::pheatmap(as.matrix(devToPlot_CM),
                              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                              border_color = NA,
                              #color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
                              scale = "row",
                              cluster_cols = FALSE,
                              cluster_rows = TRUE,
                              show_colnames = TRUE,
                              show_rownames = TRUE,
                              annotation_col = annotCol_forHeatmap,
                              annotation_colors = annotCol_forHeatmap_colors,  
                              fontsize_row = 9,
                              fontsize = 9,
                              main = paste0("Heatmap of signif. CM line ATAC-seq Drug treatment control vs 24h vs 48h. ", sbst, "- used"))

ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_CM_line_drug_treatment_control_vs_24h_48h_Motifs_", sbst, "- used", ".png")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_CM_line_drug_treatment_control_vs_24h_48h_Motifs_", sbst, "- used", ".eps")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")

# separate SH
devToPlot_SH <- devTop %>%
  dplyr::select(long_name, 
                SH_control_ATAC_S165166_REP1, SH_control_ATAC_S165166_REP2,
                SH_24h_ATAC_S165164_REP1, SH_24h_ATAC_S165164_REP2,
                SH_48h_ATAC_S165165_REP1, SH_48h_ATAC_S165165_REP2
  ) %>%
  tibble::column_to_rownames(var = "long_name")

heatmap <- pheatmap::pheatmap(as.matrix(devToPlot_SH),
                              color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                              border_color = NA,
                              #color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
                              scale = "row",
                              cluster_cols = FALSE,
                              cluster_rows = TRUE,
                              show_colnames = TRUE,
                              show_rownames = TRUE,
                              annotation_col = annotCol_forHeatmap,
                              annotation_colors = annotCol_forHeatmap_colors,  
                              fontsize_row = 9,
                              fontsize = 9,
                              main = paste0("Heatmap of signif. SH line ATAC-seq Drug treatment control vs 24h vs 48h. ", sbst, "- used"))

ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_SH_line_drug_treatment_control_vs_24h_48h_Motifs_", sbst, "- used", ".png")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, paste0("heatmap_ATAC_SH_line_drug_treatment_control_vs_24h_48h_Motifs_", sbst, "- used", ".eps")), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")


##############################################################################################################
# combining peakIds from both control results to create a universal peak space
drg_trtm_ATAC
ATAC_dds_results_control <- c(
  ATAC_dds_results_24h$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(PeakId),
  ATAC_dds_results_48h$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% dplyr::pull(PeakId)
) %>% 
  unique()
h24_specific_peaks <- ATAC_dds_results_24h$results_signif %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(PeakId)
h48_specific_peaks <- ATAC_dds_results_48h$results_signif %>% dplyr::filter(log2FoldChange > 0) %>% dplyr::pull(PeakId)

opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["species"]] <- "9606"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- TFBSTools::getMatrixSet(JASPAR2022, opts)
TEAD_AP1_motifsToScan_names <- c("MA0090.3", "MA0808.1", "MA0809.2", "MA1121.1",
                                 "MA0099.3")
TEAD_AP1_motifsToScan <- motifsToScan[TEAD_AP1_motifsToScan_names,]

# Using different subsets of peaks to check 
drg_trtm_ATAC
counts_consensus_filt_control <- drg_trtm_ATAC[drg_trtm_ATAC@rowRanges@elementMetadata@listData[["PeakId"]] %in% ATAC_dds_results_control, ]
counts_consensus_filt_24h <- drg_trtm_ATAC[drg_trtm_ATAC@rowRanges@elementMetadata@listData[["PeakId"]] %in% h24_specific_peaks, ]
counts_consensus_filt_48h <- drg_trtm_ATAC[drg_trtm_ATAC@rowRanges@elementMetadata@listData[["PeakId"]] %in% h48_specific_peaks, ]

# Create a set of random peaks for 24h 48h and control
# Load genome
genome_used <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
seqnames(genome_used) <- sub("^.{1,2}_", "", seqnames(genome_used))
seqnames(genome_used) <- sub("v", ".", seqnames(genome_used))
seqnames(genome_used) <- sub("^M", "MT", seqnames(genome_used))
# load mask regions 
blacklist_regions <- read.delim("~/workspace/neuroblastoma/resources/hg38-blacklist.v2.bed", header = FALSE)
blacklist_regions <- makeGRangesFromDataFrame(blacklist_regions, 
                                              seqnames.field = "V1", 
                                              start.field = "V2", 
                                              end.field = "V3")

# Create set of random Peaks for control 
initial_peak_set <- reduce(counts_consensus_filt_control@rowRanges)
random_peaks_control <- initial_peak_set
random_peaks_control@seqnames <- droplevels(random_peaks_control@seqnames)
random_peaks_control <- regioneR::randomizeRegions(random_peaks_control,
                                               allow.overlaps = FALSE,
                                               genome = genome_used,
                                               per.chromosome = TRUE,
                                               mask = blacklist_regions)
random_peaks_control <- as.data.frame(random_peaks_control, row.names = NULL, optional = FALSE)
random_peaks_control$seqnames <- droplevels(random_peaks_control$seqnames)
random_peaks_control <- makeGRangesFromDataFrame(random_peaks_control)

# Create set of random Peaks for 24h 
initial_peak_set <- reduce(counts_consensus_filt_24h@rowRanges)
random_peaks_24h <- initial_peak_set
random_peaks_24h@seqnames <- droplevels(random_peaks_24h@seqnames)
random_peaks_24h <- regioneR::randomizeRegions(random_peaks_24h,
                                               allow.overlaps = FALSE,
                                               genome = genome_used,
                                               per.chromosome = TRUE,
                                               mask = blacklist_regions)
random_peaks_24h <- as.data.frame(random_peaks_24h, row.names = NULL, optional = FALSE)
random_peaks_24h$seqnames <- droplevels(random_peaks_24h$seqnames)
random_peaks_24h <- makeGRangesFromDataFrame(random_peaks_24h)


# Create set of random Peaks for 48h 
initial_peak_set <- reduce(counts_consensus_filt_48h@rowRanges)
random_peaks_48h <- initial_peak_set
random_peaks_48h@seqnames <- droplevels(random_peaks_48h@seqnames)
random_peaks_48h <- regioneR::randomizeRegions(random_peaks_48h,
                                               allow.overlaps = FALSE,
                                               genome = genome_used,
                                               per.chromosome = TRUE,
                                               mask = blacklist_regions)
random_peaks_48h <- as.data.frame(random_peaks_48h, row.names = NULL, optional = FALSE)
random_peaks_48h$seqnames <- droplevels(random_peaks_48h$seqnames)
random_peaks_48h <- makeGRangesFromDataFrame(random_peaks_48h)



counts_consensus_filt_control <- chromVAR::addGCBias(counts_consensus_filt_control, genome = BSgenome.Hsapiens.UCSC.hg38) 
counts_consensus_filt_24h <- chromVAR::addGCBias(counts_consensus_filt_24h, genome = BSgenome.Hsapiens.UCSC.hg38) 
counts_consensus_filt_48h <- chromVAR::addGCBias(counts_consensus_filt_48h, genome = BSgenome.Hsapiens.UCSC.hg38) 

# For random peaks we don't do GC correction

# Having corrected for bias, we can use the matchMotifs function to identify motifs under our ATACseq peaks.
# Here we supply our RangedSummarizedExperiment of counts in peaks and the genome of interest to the matchMotifs function and use the default out of matches.

# find TEAD and AP1 motifs in real control data
motif_matches_control <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                              subject = counts_consensus_filt_control, 
                                              genome = BSgenome.Hsapiens.UCSC.hg38, 
                                              out = "positions")
control_TEAD <- c(motif_matches_control$MA0090.3, motif_matches_control$MA0808.1, motif_matches_control$MA0809.2, motif_matches_control$MA1121.1)
control_TEAD <- reduce(control_TEAD)
control_AP1 <- reduce(motif_matches_control$MA0099.3)
control_AP1_TEAD_overlaps <- findOverlaps(control_AP1, control_TEAD, ignore.strand = TRUE, maxgap = 50)

# find TEAD and AP1 motifs in simulated MES data
motif_matches_control_random <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                                     subject = random_peaks_control, 
                                                     genome = BSgenome.Hsapiens.UCSC.hg38, 
                                                     out = "positions")
control_random_TEAD <- c(motif_matches_control_random$MA0090.3, motif_matches_control_random$MA0808.1, motif_matches_control_random$MA0809.2, motif_matches_control_random$MA1121.1)
control_random_TEAD <- reduce(control_random_TEAD)
control_random_AP1 <- reduce(motif_matches_control_random$MA0099.3)
control_random_AP1_TEAD_overlaps <- findOverlaps(control_random_AP1, control_random_TEAD, ignore.strand = TRUE, maxgap = 50)

# find TEAD and AP1 motifs in real 24h data
motif_matches_24h <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                              subject = counts_consensus_filt_24h, 
                                              genome = BSgenome.Hsapiens.UCSC.hg38, 
                                              out = "positions")
h24_TEAD <- c(motif_matches_24h$MA0090.3, motif_matches_24h$MA0808.1, motif_matches_24h$MA0809.2, motif_matches_24h$MA1121.1)
h24_TEAD <- reduce(h24_TEAD)
h24_AP1 <- reduce(motif_matches_24h$MA0099.3)
h24_AP1_TEAD_overlaps <- findOverlaps(h24_AP1, h24_TEAD, ignore.strand = TRUE, maxgap = 50)

# find TEAD and AP1 motifs in simulated 24h data
motif_matches_24h_random <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                                     subject = random_peaks_24h, 
                                                     genome = BSgenome.Hsapiens.UCSC.hg38, 
                                                     out = "positions")
h24_random_TEAD <- c(motif_matches_24h_random$MA0090.3, motif_matches_24h_random$MA0808.1, motif_matches_24h_random$MA0809.2, motif_matches_24h_random$MA1121.1)
h24_random_TEAD <- reduce(h24_random_TEAD)
h24_random_AP1 <- reduce(motif_matches_24h_random$MA0099.3)
h24_random_AP1_TEAD_overlaps <- findOverlaps(h24_random_AP1, h24_random_TEAD, ignore.strand = TRUE, maxgap = 50)


# find TEAD and AP1 motifs in real 48h data
motif_matches_48h <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                              subject = counts_consensus_filt_48h, 
                                              genome = BSgenome.Hsapiens.UCSC.hg38, 
                                              out = "positions")
h48_TEAD <- c(motif_matches_48h$MA0090.3, motif_matches_48h$MA0808.1, motif_matches_48h$MA0809.2, motif_matches_48h$MA1121.1)
h48_TEAD <- reduce(h48_TEAD)
h48_AP1 <- reduce(motif_matches_48h$MA0099.3)
h48_AP1_TEAD_overlaps <- findOverlaps(h48_AP1, h48_TEAD, ignore.strand = TRUE, maxgap = 50)

# find TEAD and AP1 motifs in simulated 48h data
motif_matches_48h_random <- motifmatchr::matchMotifs(pwms = TEAD_AP1_motifsToScan, 
                                                     subject = random_peaks_48h, 
                                                     genome = BSgenome.Hsapiens.UCSC.hg38, 
                                                     out = "positions")
h48_random_TEAD <- c(motif_matches_48h_random$MA0090.3, motif_matches_48h_random$MA0808.1, motif_matches_48h_random$MA0809.2, motif_matches_48h_random$MA1121.1)
h48_random_TEAD <- reduce(h48_random_TEAD)
h48_random_AP1 <- reduce(motif_matches_48h_random$MA0099.3)
h48_random_AP1_TEAD_overlaps <- findOverlaps(h48_random_AP1, h48_random_TEAD, ignore.strand = TRUE, maxgap = 50)


summary_table_peaks <- data.frame(total_number_of_peaks = c(length(counts_consensus_filt_control), 
                                                            length(random_peaks_control),
                                                            length(counts_consensus_filt_24h),
                                                            length(random_peaks_24h),
                                                            length(counts_consensus_filt_48h),
                                                            length(random_peaks_48h)
                                                            ),
                                  TEAD_sites = c(length(control_TEAD),
                                                 length(control_random_TEAD),
                                                 length(h24_TEAD),
                                                 length(h24_random_TEAD),
                                                 length(h48_TEAD),
                                                 length(h48_random_TEAD)
                                                 ),
                                  AP1_sites = c(length(control_AP1),
                                                length(control_random_AP1),
                                                length(h24_AP1),
                                                length(h24_random_AP1),
                                                length(h48_AP1),
                                                length(h48_random_AP1)
                                                ),
                                  TEAD_AP1_coloc_sites = c(length(control_AP1_TEAD_overlaps),
                                                           length(control_random_AP1_TEAD_overlaps),
                                                           length(h24_AP1_TEAD_overlaps),
                                                           length(h24_random_AP1_TEAD_overlaps),
                                                           length(h48_AP1_TEAD_overlaps),
                                                           length(h48_random_AP1_TEAD_overlaps)
                                                           ),
                                  row.names = c("control", "Simulated control", "h24", "Simulated h24", "h48", "Simulated h48")
)
summary_table_peaks


# Create contingency tables to do Ftest on overlap 24h control
contingency_table <- matrix(c(length(h24_AP1_TEAD_overlaps), length(h24_TEAD) + length(h24_AP1) - length(h24_AP1_TEAD_overlaps), 
                              length(control_AP1_TEAD_overlaps), length(control_TEAD) + length(control_AP1) - length(control_AP1_TEAD_overlaps)
), nrow=2, byrow=TRUE)
fisher_results <- fisher.test(contingency_table, alternative="two.sided")
fisher_summary_table <- data.table::data.table(comparison = "MES_vs_ADR",
                                               pval = fisher_results$p.value,
                                               odds_ratio = fisher_results$estimate)

# Fisher test comparing real MES AP1_TEAD peaks vs simulated
contingency_table_MES_vs_simulated <- matrix(c(length(MES_AP1_TEAD_overlaps), length(MES_TEAD) + length(MES_AP1) - length(MES_AP1_TEAD_overlaps),
                                               length(MES_random_AP1_TEAD_overlaps), length(MES_random_TEAD) + length(MES_random_AP1) - length(MES_random_AP1_TEAD_overlaps)
), nrow=2, byrow=TRUE)
row.names(contingency_table_MES_vs_simulated) <- c("MES_peaks", "MES_peaks_random")
colnames(contingency_table_MES_vs_simulated) <- c("TEAD_AB1_colocalisation_present", "TEAD_AB1_colocalisation_absent")
contingency_table_MES_vs_simulated
fisher_results <- fisher.test(contingency_table_MES_vs_simulated, alternative="two.sided")
fisher_summary_table <- rbind(fisher_summary_table,
                              data.table::data.table(comparison = "MES_vs_MESsim",
                                                     pval = fisher_results$p.value,
                                                     odds_ratio = fisher_results$estimate))

