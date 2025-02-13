###
# Title: RNA-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
### 

# importing only key functions that are actually used - not to polute namespace!
import::from(.from = DESeq2, .all=TRUE)
import::from(openxlsx, createWorkbook, addWorksheet, writeData, saveWorkbook)
import::from(.from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"), 
             "filterDatasets",
             "generatePCA", 
             "generateEnsemblAnnotation", 
             "generateResults",
             "meanExprsPerGroup",
             "plotVolcano",
             "generatePCA_repel",
             .character_only=TRUE) # used for filtering
library(ggplot2)
library(GSVA)
library(pheatmap)
library(xcore)
library(ExperimentHub)
library(xcoredata)
library(stringr)
library(org.Hs.eg.db)
library(fgsea)
library(dplyr)

# load annotation table. Could use Biomart, alternatively, but had this 
param_list <- list(
  abs_filt_samples=2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 0.58,
  var_expl_needed = 0.6,
  biomart_host="http://www.ensembl.org", 
  biomart_dataset="hsapiens_gene_ensembl", 
  biomart_Ens_version="Ensembl Genes 109"
)

# load gene sets
gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean = TRUE)
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean = TRUE)
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean = TRUE)
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean = TRUE)
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean = TRUE)
gs_C6_onco <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C6"), clean = TRUE)
gs_wang_hippo <- list(wang_hippo_set = c("CCN1", "CCN2", "AMOTL2", "ANKRD1", "IGFBP3", "F3", "FJX1", "NUAK2", "LATS2", "CRIM1", "GADD45A",
                                         "TGFB2", "PTPN14", "NT5E", "FOXF2", "AXL", "DOCK5", "ASAP1", "RBMS3", "MYOF", "ARHGEF17", "CCDC80"))

RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)



# load annotation
path_folder_rna_seq_data <- "~/workspace/neuroblastoma/data/RNAseq"
annotationData <- read.table(
  file = file.path(path_folder_rna_seq_data, "rnaseq_deseq_global_annotation_gene.tsv"),
  sep = "\t",
  header = TRUE
)
# This is necessary to rename columns so that extraction of the data works correctly
colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")



deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_yap_taz_inhibition/"

#loading the new DESeq2 file
load("~/workspace/neuroblastoma/data/RNAseq/taz_yap_inhibition/deseq2.dds.RData")
dds$cell_line <- factor(dds$Group1, levels = c("CM", "SH"))
dds$timepoint <- factor(dds$Group2, levels = c("control", "24h", "48h"))
dds$replicate <- dds$Group3

design(dds) <- as.formula("~ cell_line + timepoint")

dds <- filterDatasets(dds, 
                      abs_filt = TRUE, 
                      abs_filt_samples = param_list$abs_filt_samples)
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds)

#stabilize variance
vsd <- DESeq2::vst(dds, blind = TRUE) # blind = TRUE for QC
#generate PCA plots
pca_deg <- generatePCA_repel(transf_object = vsd, 
                             cond_interest_varPart = c("cell_line","timepoint"), 
                             color_variable = "timepoint", 
                             shape_variable = "cell_line",
                             ntop_genes = 1000) +
  ggtitle("Original dataset") 
pdf(file = file.path(deg_dir, "pca_original_dataset.pdf"))
pca_deg
dev.off()


# Batch correction
transf_batch_NObatch_experiment <- vsd
transf_batch_NObatch_experiment_count <- limma::removeBatchEffect(SummarizedExperiment::assay(transf_batch_NObatch_experiment),
                                                                  transf_batch_NObatch_experiment$cell_line)
SummarizedExperiment::assay(transf_batch_NObatch_experiment) <- transf_batch_NObatch_experiment_count
pca_deg_NObatch_experiment <- generatePCA_repel(transf_object = transf_batch_NObatch_experiment, 
                                                cond_interest_varPart = c("cell_line","timepoint"), 
                                                color_variable = "timepoint", 
                                                shape_variable = "cell_line",
                                                ntop_genes = 1000) +
  ggtitle("Original dataset - YAP1/TAZ1 inhibition; cell_line batch correction") 

pdf(file = file.path(deg_dir, "pca_with_corrected_batch_effect_from_cell_line.pdf"))
pca_deg_NObatch_experiment
dev.off()

resultsNames(dds)


######################################
## 24H vs control
resultsNames(dds)
deg_results <- generateResults(
  dds_object = dds,
  coeff_name = "timepoint_24h_vs_control",
  cond_numerator = "24h",
  cond_denominator = "control",
  cond_variable = "timepoint",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

XLSX_OUT <- createWorkbook()
addWorksheet(XLSX_OUT, "results_signif")
addWorksheet(XLSX_OUT, "de_details")
addWorksheet(XLSX_OUT, "results_all")

writeData(XLSX_OUT, x = deg_results$results_signif, sheet = "results_signif")
writeData(XLSX_OUT, x = deg_results$de_details, sheet = "de_details")
writeData(XLSX_OUT, x = deg_results$results_all, sheet = "results_all")

saveWorkbook(XLSX_OUT, 
             file.path(deg_dir, "cell_type_24р_vs_control.xlsx"),
             overwrite = TRUE)

# VOLCANO PLOT 24h
genes_to_highlight <- c(
  "VIM",
  "YAP1",
  "WWTR1",
  "JUN",
  "FOSL1",
  "FOSL2",
  "PHOX2B",
  "HAND2",
  "GATA3"
)
volcano_plot <- plotVolcano(
  dds_results_obj = deg_results$results_all,
  genes_of_interest = genes_to_highlight,
  plot_title = "24H vs control"
)
ggsave(
  filename = paste0(deg_dir, "24H_vs_control_volcano_plot.png"),
  plot = volcano_plot,
  width = 20, height = 20, units = "cm"
)


# Prepare data for heatmaps 24H
metadata_heatmap <- as.data.frame(colData(dds))
dds_signif <- deg_results$results_signif

heatmap_counts <- SummarizedExperiment::assay(vsd)
heatmap_counts <- SummarizedExperiment::assay(transf_batch_NObatch_experiment_count)
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% dds_signif$ensembl_id, ]

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, timepoint, sample) %>% 
  dplyr::arrange(timepoint, cell_line)
heatmap_counts <- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]

# Alternative color schemes
# color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012", "#9B2226")
# color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
# color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
color.scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
ann_colors = list(
  cell_line = c(CM = "#005f73", SH = "#FF9F1C"),
  timepoint = c(control = "#E9D8A6", `24h` = "#D9BE6D", `48h` = "#d6ac2f")
)

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "24H vs control samples. VST batch corrected for cell_line",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10)

# Produce GSEA plots
dwn_deg_results <- deg_results$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% pull(gene_symbol)

C5_GOBP <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOBP,
  test = "hypergeometric",
  background = nrow(dds)
)
C5_GOCC <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOCC,
  test = "hypergeometric",
  background = nrow(dds)
)
C5_GOMF <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOMF,
  test = "hypergeometric",
  background = nrow(dds)
)
C2_kegg <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C2_kegg,
  test = "hypergeometric",
  background = nrow(dds)
)
C2_reactome <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C2_reactome,
  test = "hypergeometric",
  background = nrow(dds)
)
C6_onco <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C6_onco,
  test = "hypergeometric",
  background = nrow(dds)
)
wang_hippo <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_wang_hippo,
  test = "hypergeometric",
  background = nrow(dds)
)

C5_GOBP_plot <- hypeR::hyp_dots(C5_GOBP, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOBP: 24h vs control downregulated") + theme_bw()
C5_GOCC_plot <- hypeR::hyp_dots(C5_GOCC, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOCC: 24h vs control downregulated") + theme_bw()
C5_GOMF_plot <- hypeR::hyp_dots(C5_GOMF, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOMF: 24h vs control downregulated") + theme_bw()
C2_kegg_plot <- hypeR::hyp_dots(C2_kegg, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "KEGG: 24h vs control downregulated") + theme_bw()
C2_rctm_plot <- hypeR::hyp_dots(C2_reactome, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "REACTOME: 24h vs control downregulated") + theme_bw()
C6_onco_plot <- hypeR::hyp_dots(C6_onco, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "ONCO: 24h vs control downregulated") + theme_bw()




######################################
## 48H vs control
resultsNames(dds)
deg_results <- generateResults(
  dds_object = dds,
  coeff_name = "timepoint_48h_vs_control",
  cond_numerator = "48h",
  cond_denominator = "control",
  cond_variable = "timepoint",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

XLSX_OUT <- createWorkbook()
addWorksheet(XLSX_OUT, "results_signif")
addWorksheet(XLSX_OUT, "de_details")
addWorksheet(XLSX_OUT, "results_all")

writeData(XLSX_OUT, x = deg_results$results_signif, sheet = "results_signif")
writeData(XLSX_OUT, x = deg_results$de_details, sheet = "de_details")
writeData(XLSX_OUT, x = deg_results$results_all, sheet = "results_all")

saveWorkbook(XLSX_OUT, 
             file.path(deg_dir, "cell_type_48h_vs_control.xlsx"),
             overwrite = TRUE)

# VOLCANO PLOT 48h
genes_to_highlight <- c(
  "VIM",
  "YAP1",
  "WWTR1",
  "JUN",
  "FOSL1",
  "FOSL2",
  "PHOX2B",
  "HAND2",
  "GATA3"
)
volcano_plot <- plotVolcano(
  dds_results_obj = deg_results$results_all,
  genes_of_interest = genes_to_highlight,
  plot_title = "48H vs control"
)
ggsave(
  filename = paste0(deg_dir, "48H_vs_control_volcano_plot.png"),
  plot = volcano_plot,
  width = 20, height = 20, units = "cm"
)


# Prepare data for heatmaps 48H
metadata_heatmap <- as.data.frame(colData(dds))
dds_signif <- deg_results$results_signif

heatmap_counts <- SummarizedExperiment::assay(vsd)
heatmap_counts <- SummarizedExperiment::assay(transf_batch_NObatch_experiment_count)
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% dds_signif$ensembl_id, ]

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, timepoint, sample) %>% 
  dplyr::arrange(timepoint, cell_line)
heatmap_counts <- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]

# Alternative color schemes
# color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012", "#9B2226")
# color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
# color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
color.scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
ann_colors = list(
  cell_line = c(CM = "#005f73", SH = "#FF9F1C"),
  timepoint = c(control = "#E9D8A6", `24h` = "#D9BE6D", `48h` = "#d6ac2f")
)

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "48H vs control samples. VST batch corrected for cell_line",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10)

# Produce GSEA plots
dwn_deg_results <- deg_results$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% pull(gene_symbol)

C5_GOBP <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOBP,
  test = "hypergeometric",
  background = nrow(dds)
)
C5_GOCC <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOCC,
  test = "hypergeometric",
  background = nrow(dds)
)
C5_GOMF <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C5_GOMF,
  test = "hypergeometric",
  background = nrow(dds)
)
C2_kegg <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C2_kegg,
  test = "hypergeometric",
  background = nrow(dds)
)
C2_reactome <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C2_reactome,
  test = "hypergeometric",
  background = nrow(dds)
)
C6_onco <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_C6_onco,
  test = "hypergeometric",
  background = nrow(dds)
)
wang_hippo <- hypeR::hypeR(
  signature = dwn_deg_results,
  genesets = gs_wang_hippo,
  test = "hypergeometric",
  background = nrow(dds)
)

C5_GOBP_plot <- hypeR::hyp_dots(C5_GOBP, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOBP: 48h vs control downregulated") + theme_bw()
C5_GOCC_plot <- hypeR::hyp_dots(C5_GOCC, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOCC: 48h vs control downregulated") + theme_bw()
C5_GOMF_plot <- hypeR::hyp_dots(C5_GOMF, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOMF: 48h vs control downregulated") + theme_bw()
C2_kegg_plot <- hypeR::hyp_dots(C2_kegg, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "KEGG: 48h vs control downregulated") + theme_bw()
C2_rctm_plot <- hypeR::hyp_dots(C2_reactome, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "REACTOME: 48h vs control downregulated") + theme_bw()
C6_onco_plot <- hypeR::hyp_dots(C6_onco, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "ONCO: 48h vs control downregulated") + theme_bw()

C5_GOBP_plot
C5_GOCC_plot
C5_GOMF_plot
C2_kegg_plot
C2_rctm_plot
C6_onco_plot




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### ffGSEA to show a shift from MES to ADR identity for all samples ####
vsd_counts_matrix <- assay(vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, timepoint, sample) %>% 
  dplyr::arrange(timepoint, cell_line)

ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines[, match(rownames(annotation_col), colnames(ssgsea_mes_adr_cellines))]

ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        annotation_col = annotation_col,
                                                        annotation_colors = ann_colors,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = FALSE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)
pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap
dev.off()
