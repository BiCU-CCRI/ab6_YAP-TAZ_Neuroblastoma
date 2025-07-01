###
# Title: RNA-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
### 

# importing only key functions that are actually used - not to polute namespace!
library(DESeq2)
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
#necessary for xcoredata
if (!dir.exists("/home/rstudio/.cache/R/ExperimentHub")){
  dir.create("/home/rstudio/.cache/R/ExperimentHub")}
library(xcoredata)
library(stringr)
library(org.Hs.eg.db)
library(fgsea)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rlang)

# Declare fGSEA function 
plot_fGSEA <- function(deg_results, gene_set_list, title_prefix, save_dir, maxSize = 500, figure_extention = ".png"){
  de_genes_ranked <- deg_results
  de_genes_ranked <- de_genes_ranked %>%
    group_by(gene_symbol) %>%
    slice_max(baseMean, n = 1, with_ties = FALSE) %>%
    ungroup()
  de_genes_ranked <- de_genes_ranked %>% arrange(desc(log2FoldChange))
  de_genes_ranked <- setNames(c(de_genes_ranked$log2FoldChange), c(de_genes_ranked$gene_symbol))
 
  fgseaRes <- fgsea(
    pathways = gene_set_list,
    stats = de_genes_ranked,
    minSize = 15,
    maxSize = maxSize
  )
   
  for(gene_set_to_plot_name in names(gene_set_list)){
    # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
    gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
    
    title_values <- fgseaRes %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
    p <- plotEnrichment(pathway = gene_set_to_plot, 
                        stats = de_genes_ranked) + 
      labs(title = paste(title_prefix, "\n",
                         title_values$pathway, "\n",
                         "NES =", title_values$NES,
                         "padj = ", title_values$padj))
    plot(p)
    ggsave(file.path(save_dir, paste0(title_prefix, "_", title_values$pathway, figure_extention)), plot = p)
  }
  
}

produce_GSEA_plots <- function(gene_signature, dds_object, additional_title, path_to_pdf_report){
  C5_GOBP <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C5_GOBP,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  C5_GOCC <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C5_GOCC,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  C5_GOMF <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C5_GOMF,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  C2_kegg <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C2_kegg,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  C2_reactome <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C2_reactome,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  C6_onco <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_C6_onco,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  wang_hippo <- hypeR::hypeR(
    signature = gene_signature,
    genesets = gs_wang_hippo,
    test = "hypergeometric",
    background = nrow(dds_object)
  )
  
  C5_GOBP_plot <- hypeR::hyp_dots(C5_GOBP, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOBP: ", additional_title)) + theme_bw()
  C5_GOCC_plot <- hypeR::hyp_dots(C5_GOCC, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOCC: ", additional_title))+ theme_bw()
  C5_GOMF_plot <- hypeR::hyp_dots(C5_GOMF, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOMF: ", additional_title)) + theme_bw()
  C2_kegg_plot <- hypeR::hyp_dots(C2_kegg, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("KEGG: ", additional_title)) + theme_bw()
  C2_rctm_plot <- hypeR::hyp_dots(C2_reactome, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("REACTOME: ", additional_title)) + theme_bw()
  C6_onco_plot <- hypeR::hyp_dots(C6_onco, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("ONCO: ", additional_title)) + theme_bw()
  
  pdf(file = file.path(path_to_pdf_report, paste0(additional_title, ".pdf")))
  print(C5_GOBP_plot)
  print(C5_GOCC_plot)
  print(C5_GOMF_plot)
  print(C2_kegg_plot)
  print(C2_rctm_plot)
  print(C6_onco_plot)
  dev.off()
}

# load annotation table. Could use Biomart, alternatively, but had this 
param_list <- list(
  abs_filt_samples=2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 0.37, #MAybe need to be changed later 
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
gs_c2_PID <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:PID", clean = TRUE)

# RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
# RNA_SEQ_data <- RNA_SEQ_data %>%
#   mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES")) %>%
#   select(gene_symbol, Term)
# ADRN_signature_list <- RNA_SEQ_data %>% filter(Term == "ADRN") %>% pull(gene_symbol) %>% as.vector()
# MES_signature_list <- RNA_SEQ_data %>% filter(Term == "MES") %>% pull(gene_symbol) %>% as.vector()
AP1_maaynlab <- read.csv("~/workspace/neuroblastoma/resources/ap1_gene_list", header = FALSE) %>% pull(V1) %>% as.vector()

gene_set_list <- list(
  GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
  REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
  C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
  Hippo_Wang = gs_wang_hippo$wang_hippo_set,
  C2_PID_pathway = gs_c2_PID[["genesets"]][["Ap1 Pathway"]],
  AP1_maaynlab = AP1_maaynlab
)

RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)
# Adding strict filter
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 3)
mes_adrn_gene_list_strict <- RNA_SEQ_data %>%
  filter(padj < param_list$padj_cutoff) %>% 
  mutate(Term = case_when(
    FoldChange < -10 ~ "ADRN",
    FoldChange >  10 ~ "MES",
    TRUE ~ "no_term"
))
sig_list_our_data_strict <- list(
  Strict_Aderenergic = mes_adrn_gene_list_strict$gene_symbol[mes_adrn_gene_list_strict$Term == "ADRN"],
  Strict_Mesenchymal = mes_adrn_gene_list_strict$gene_symbol[mes_adrn_gene_list_strict$Term == "MES"]
)

# loading signatures from Groeningen paper
signatures_raw <- readr::read_tsv(file = "/home/rstudio/workspace/neuroblastoma/resources/mes_adr_ncc_noradr_Groeningen_genesets.tsv",
                                  skip_empty_rows = TRUE) # do not add na for empty
sig_list_groeningen <- list(
  #NCC_like = signatures_raw$`NCC-like`[!is.na(signatures_raw$`NCC-like`)],
  #Noradrenergic = signatures_raw$Noradrenergic[!is.na(signatures_raw$Noradrenergic)],
  Groen_Mesenchymal  = signatures_raw$Mesenchymal[!is.na(signatures_raw$Mesenchymal)],
  Groen_Adrenergic = signatures_raw$Adrenergic[!is.na(signatures_raw$Adrenergic)]
)
lapply(sig_list_groeningen, length)

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
if(!dir.exists(deg_dir)){
  dir.create(deg_dir)
} 

#loading the new DESeq2 file
load("~/workspace/neuroblastoma/data/RNAseq/taz_yap_inhibition/deseq2.dds.RData")
dds$cell_line <- factor(dds$Group1, levels = c("CM", "SH"))
dds$timepoint <- factor(dds$Group2, levels = c("control", "24h", "48h"))
dds$replicate <- dds$Group3
comparison_groups <- paste0(dds$cell_line, "_", dds$timepoint)
dds$comparison_group <- factor(comparison_groups, levels = unique(comparison_groups))


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

# Barplots and statistics for CTGF ANKRD1
genes_of_interest_symbols <- c("CCN2", "ANKRD1")
annotationData$ensembl_id[annotationData$gene_symbol %in% genes_of_interest_symbols]
genes_of_interest <- annotationData$ensembl_id[annotationData$gene_symbol %in% c("CCN2", "ANKRD1")]
dds_sub <- dds[rownames(dds) %in% genes_of_interest, ]
norm_counts <- counts(dds_sub, normalized = TRUE) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

norm_counts <- left_join(norm_counts, as.data.frame(colData(dds)) , by = "sample")

norm_counts$timepoint <- factor(norm_counts$timepoint, levels = c("control", "24h", "48h"))
norm_counts$cell_line <- factor(norm_counts$cell_line)
norm_counts$gene <- factor(norm_counts$gene)
norm_counts <- norm_counts %>% mutate(gene = recode(gene, 
                                                    "ENSG00000118523" = "CCN2",
                                                    "ENSG00000148677" = "ANKRD1"))
# Generate plot
for (gene in genes_of_interest_symbols){
  for (cell_line in c("CM", "SH")) {
   p <- norm_counts %>% dplyr::filter(gene == !!gene, cell_line == !!cell_line) %>%
      ggplot(aes(x = timepoint, y = expression, fill = timepoint)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      stat_compare_means(method = "t.test", 
                         comparisons = list(c("control", "24h"), c("control", "48h")),
                        label = "p.format") +
      theme_minimal(base_size = 14) +
      labs(
        title = paste0("Gene Expression by Timepoint in ", cell_line, " cell line"),
        x = "Timepoint",
        y = paste0("Normalized expression \n", gene )
      ) +
      scale_fill_brewer(palette = "Set2")
   plot(p)
  }
}


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
             file.path(deg_dir, "cell_type_24h_vs_control.xlsx"),
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
  log2FC_cutoff = param_list$log2FC_cutoff,
  padj_cutoff = param_list$padj_cutoff,
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
produce_GSEA_plots(gene_signature = dwn_deg_results, 
                   dds_object = dds, 
                   additional_title = "24h vs control downregulated", 
                   path_to_pdf_report = deg_dir)

# deg_results$results_all <- deg_results$results_all %>% filter(gene_biotype == "protein_coding")
# deg_results$results_all <- deg_results$results_all %>% filter(is.finite(log2FoldChange))

# Produce fGSEA
plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = gene_set_list,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data_strict,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_groeningen,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)
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
  log2FC_cutoff = param_list$log2FC_cutoff,
  padj_cutoff = param_list$padj_cutoff,
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
produce_GSEA_plots(gene_signature = dwn_deg_results, 
                   dds_object = dds, 
                   additional_title = "48h vs control downregulated", 
                   path_to_pdf_report = deg_dir)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = gene_set_list,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
           )

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data_strict,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_groeningen,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### GSVA to show a shift from MES to ADR identity for all samples ####
vsd_counts_matrix <- assay(vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)

metadata_heatmap <- as.data.frame(colData(dds))
annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, timepoint, sample) %>% 
  dplyr::arrange(timepoint, cell_line)
ann_colors = list(
  cell_line = c(CM = "#005f73", SH = "#FF9F1C"),
  timepoint = c(control = "#E9D8A6", `24h` = "#D9BE6D", `48h` = "#d6ac2f")
)

ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines[, match(rownames(annotation_col), colnames(ssgsea_mes_adr_cellines))]

ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        annotation_col = annotation_col,
                                                        annotation_colors = ann_colors,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = FALSE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)

ggsave(
  filename = paste0(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"),
  plot = ssgsea_mes_adr_ncc_noradr_heatmap,
  width = 20, height = 10, units = "cm"
)



# MES separately
ssgsea_mes_adr_cellines_MES <- as.data.frame(ssgsea_mes_adr_cellines)
ssgsea_mes_adr_cellines_MES <- ssgsea_mes_adr_cellines_MES["Mesenchymal",]
ssgsea_mes_adr_ncc_noradr_heatmap_MES <- pheatmap::pheatmap(ssgsea_mes_adr_cellines_MES,
                                                        scale = "row",
                                                        annotation_col = annotation_col,
                                                        annotation_colors = ann_colors,
                                                        cluster_rows = FALSE,
                                                        cluster_cols = FALSE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)
ggsave(
  filename = paste0(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ_MES.pdf"),
  plot = ssgsea_mes_adr_ncc_noradr_heatmap_MES,
  width = 20, height = 10, units = "cm"
)


# ADR separately
ssgsea_mes_adr_cellines_ADR <- as.data.frame(ssgsea_mes_adr_cellines)
ssgsea_mes_adr_cellines_ADR <- ssgsea_mes_adr_cellines_ADR["Aderenergic",]
ssgsea_mes_adr_ncc_noradr_heatmap_ADR <- pheatmap::pheatmap(ssgsea_mes_adr_cellines_ADR,
                                                            scale = "row",
                                                            annotation_col = annotation_col,
                                                            annotation_colors = ann_colors,
                                                            cluster_rows = FALSE,
                                                            cluster_cols = FALSE,
                                                            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                            show_colnames = TRUE)
ggsave(
  filename = paste0(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ_ADR.pdf"),
  plot = ssgsea_mes_adr_ncc_noradr_heatmap_ADR,
  width = 20, height = 10, units = "cm"
)

# # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # #
# now do comparisons for CM and SH separately
# also separating 24h vs control and 48h vs control
load("~/workspace/neuroblastoma/data/RNAseq/taz_yap_inhibition/deseq2.dds.RData")
dds$cell_line <- factor(dds$Group1, levels = c("CM", "SH"))
dds$timepoint <- factor(dds$Group2, levels = c("control", "24h", "48h"))
dds$replicate <- dds$Group3
comparison_groups <- paste0(dds$cell_line, "_", dds$timepoint)
dds$comparison_group <- factor(comparison_groups, 
                               levels = c("CM_control", "CM_24h", "CM_48h", 
                                          "SH_control", "SH_24h", "SH_48h" ))


for (cell_line_name in c("CM", "SH")){
  #cell_line_name <- "CM"
  
  print(paste0("Processing: ", cell_line_name))
  cols_to_subset <- startsWith(colnames(dds), cell_line_name)
  dds_cell_line_subset <- dds[, cols_to_subset]
  dds_cell_line_subset$comparison_group <- droplevels(dds_cell_line_subset$comparison_group)
  
  design(dds_cell_line_subset) <- as.formula("~ comparison_group")
  
  dds_cell_line_subset <- filterDatasets(dds_cell_line_subset, 
                                         abs_filt = TRUE, 
                                         abs_filt_samples = param_list$abs_filt_samples)
  dds_cell_line_subset <- DESeq2::estimateSizeFactors(dds_cell_line_subset)
  dds_cell_line_subset <- DESeq2::DESeq(dds_cell_line_subset)
  
  for(comparison_name in resultsNames(dds_cell_line_subset)[2:3]){
    #comparison_name <- resultsNames(dds_cell_line_subset)[2]
    print(paste0("Processing: ", comparison_name))
    group_parameters <- str_split(comparison_name, "_", simplify = TRUE)[c(3,4,7)]
    cond_numerator <- paste0(group_parameters[1], "_", group_parameters[2])
    cond_denominator <- paste0(group_parameters[1], "_", group_parameters[3])
    
    deg_results <- generateResults(
      dds_object = dds_cell_line_subset,  
      coeff_name = comparison_name,
      cond_numerator = cond_numerator,
      cond_denominator = cond_denominator,
      cond_variable = "comparison_group",
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
                 file.path(deg_dir, paste0(comparison_name, "_", ".xlsx")),
                 overwrite = TRUE)
    
    # VOLCANO PLOT 
    genes_to_highlight <- c(
      "VIM",
      "YAP1",
      "WWTR1",
      "JUN",
      "FOSL1",
      "FOSL2",
      "PHOX2B",
      "HAND2",
      "GATA3",
      "NOTCH1"
    )
    volcano_plot <- plotVolcano(
      dds_results_obj = deg_results$results_all,
      genes_of_interest = genes_to_highlight,
      plot_title = comparison_name, 
      log2FC_cutoff = 0.37, 
      padj_cutoff = 0.05
    ) + 
      scale_x_continuous(limits = c(-5, 7.5)) +
      scale_y_continuous(limits = c(-5, 320))
    ggsave(
      filename = paste0(deg_dir, paste0(comparison_name, "_more_highlited_genes.pdf")),
      plot = volcano_plot,
      width = 20, height = 20, units = "cm"
    )
    
    # VOLCANO PLOT with less genes
    genes_to_highlight <- c(
      "FOSL1",
      "NOTCH1"
    )
    volcano_plot <- plotVolcano(
      dds_results_obj = deg_results$results_all,
      genes_of_interest = genes_to_highlight,
      plot_title = comparison_name, 
      log2FC_cutoff = 0.37, 
      padj_cutoff = 0.05
    ) + 
      scale_x_continuous(limits = c(-5, 7.5)) +
      scale_y_continuous(limits = c(-5, 320))
    ggsave(
      filename = paste0(deg_dir, paste0(comparison_name, "_less_highlited_genes.pdf")),
      plot = volcano_plot,
      width = 20, height = 20, units = "cm"
    )
  }
 
}







# # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # 
# Analysis of JunDN #####
deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_yap_taz_inhibition/Jun_DN/"
if(!dir.exists(deg_dir)){
  dir.create(deg_dir)
} 

unzip(zipfile = "~/workspace/neuroblastoma/data/RNAseq/JunDN/rnaseq_deseq_jundn_counts_raw.tsv.zip",
      overwrite = TRUE,
      exdir = "~/workspace/neuroblastoma/data/RNAseq/JunDN/")

raw_JunDN_counts <- read.csv(file = "~/workspace/neuroblastoma/data/RNAseq/JunDN/rnaseq_deseq_jundn_counts_raw.tsv", 
                             header = TRUE, 
                             sep = "\t")
raw_JunDN_counts$gene_id <- stringr::str_replace(raw_JunDN_counts$gene_id, pattern = "\\..*", replacement = "") 
row.names(raw_JunDN_counts) <- raw_JunDN_counts$gene_id
raw_JunDN_counts <- raw_JunDN_counts %>% select(starts_with("JunDN"))

# remove JUN from the count table as it's a DN experiment
raw_JunDN_counts <- raw_JunDN_counts %>% filter(row.names(.) != "ENSG00000177606")

coldata <- data.frame(
  sample = c("JunDN_ctrl_2_RNA_S165182", 
             "JunDN_ctrl_3_RNA_S165183", 
             "JunDN_dox_1_RNA_S165180", 
             "JunDN_dox_2_RNA_S165181", 
             "JunDN_dox_3_RNA_S165186"),
  replicate = c("2","3","1","2","3"),
  group = c("ctrl", "ctrl", "dox", "dox", "dox")
)
row.names(coldata) <- coldata$sample
JunDN_dds <- DESeqDataSetFromMatrix(as.matrix(raw_JunDN_counts),
                                   colData = coldata,
                                   design = as.formula(~group))

JunDN_dds <- filterDatasets(JunDN_dds, 
                      abs_filt = TRUE, 
                      abs_filt_samples = param_list$abs_filt_samples)
JunDN_dds <- DESeq2::estimateSizeFactors(JunDN_dds)
JunDN_dds <- DESeq2::DESeq(JunDN_dds)

#stabilize variance
JunDN_vsd <- DESeq2::vst(JunDN_dds, blind = TRUE) # blind = TRUE for QC
#generate PCA plots
pca_deg <- generatePCA_repel(transf_object = JunDN_vsd, 
                             cond_interest_varPart = c("group", "replicate"), 
                             color_variable = "group", 
                             shape_variable = "replicate",
                             ntop_genes = 1000) +
  ggtitle("Original dataset") 
pdf(file = file.path(deg_dir, "pca_original_dataset.pdf"))
pca_deg
dev.off()




######################################
## DOX vs Control
resultsNames(JunDN_dds)

deg_results <- generateResults(
  dds_object = JunDN_dds,
  coeff_name = "group_dox_vs_ctrl",
  cond_numerator = "dox",
  cond_denominator = "ctrl",
  cond_variable = "group",
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
             file.path(deg_dir, "JunDN_dox_vs_ctrl.xlsx"),
             overwrite = TRUE)

# VOLCANO PLOT 
genes_to_highlight <- c(
  "VIM",
  "YAP1",
  "WWTR1",
  "JUN",
  "FOSL1",
  "FOSL2",
  "PHOX2B",
  "HAND2",
  "GATA3",
  "NOTCH1"
)
volcano_plot <- plotVolcano(
  dds_results_obj = deg_results$results_all,
  genes_of_interest = genes_to_highlight,
  plot_title = "Jun DN dox vs control", 
  log2FC_cutoff = 0.37, 
  padj_cutoff = 0.05
) + 
  scale_x_continuous(limits = c(-5, 7.5)) +
  scale_y_continuous(limits = c(-5, 320))

ggsave(
  filename = paste0(deg_dir, "JunDN_dox_vs_ctrl_more_highlited_genes.pdf"),
  plot = volcano_plot,
  width = 20, height = 20, units = "cm"
)

# VOLCANO PLOT with less genes
genes_to_highlight <- c(
  "FOSL1",
  "NOTCH1"
)
volcano_plot <- plotVolcano(
  dds_results_obj = deg_results$results_all,
  genes_of_interest = genes_to_highlight,
  plot_title = "Jun DN dox vs control", 
  log2FC_cutoff = 0.37, 
  padj_cutoff = 0.05
) + 
  scale_x_continuous(limits = c(-5, 7.5)) +
  scale_y_continuous(limits = c(-5, 320))
ggsave(
  filename = paste0(deg_dir, "JunDN_dox_vs_ctrl_less_highlited_genes.pdf"),
  plot = volcano_plot,
  width = 20, height = 20, units = "cm"
)

# Prepare data for JunDN heatmap - show all DEGs
metadata_heatmap <- as.data.frame(colData(JunDN_dds))
dds_signif <- deg_results$results_signif

heatmap_counts <- SummarizedExperiment::assay(JunDN_vsd)
#heatmap_counts <- SummarizedExperiment::assay(transf_batch_NObatch_experiment_count)
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% dds_signif$ensembl_id, ]

annotation_col <- metadata_heatmap %>%
  dplyr::select(group, sample) %>% 
  dplyr::arrange(group)
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
                              main = "JunDN dox vs ctrl",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10)


# SHow only genes that belong to MES ADR signatures
metadata_heatmap <- as.data.frame(colData(JunDN_dds))
dds_signif <- deg_results$results_signif

heatmap_counts <- SummarizedExperiment::assay(JunDN_vsd)
#heatmap_counts <- SummarizedExperiment::assay(transf_batch_NObatch_experiment_count)
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% dds_signif$ensembl_id, ]
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% mes_adrn_gene_list$ensembl_id, ]
annotation_row <- mes_adrn_gene_list %>% select(ensembl_id, Term) %>% filter(ensembl_id %in% rownames(heatmap_counts))
row.names(annotation_row) <- annotation_row$ensembl_id
annotation_row <- annotation_row %>% select(Term)

annotation_col <- metadata_heatmap %>%
  dplyr::select(group, sample) %>% 
  dplyr::arrange(group)
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
                              main = "JunDN dox vs ctrl",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              annotation_row = annotation_row,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10)

# order counts by adr mes


heatmap <- pheatmap::pheatmap(heatmap_counts[row.names(annotation_row %>% arrange(Term)), ],
                              main = "JunDN dox vs ctrl",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              annotation_row = annotation_row,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              cluster_rows = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10)

# Produce fGSEA plots
dwn_deg_results <- deg_results$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% pull(gene_symbol)


produce_GSEA_plots(gene_signature = dwn_deg_results,
                   dds_object = JunDN_dds, 
                   additional_title = "JunDN_GSEA_downregulated_", 
                   path_to_pdf_report = deg_dir)

# deg_results$results_all <- deg_results$results_all %>% filter(gene_biotype == "protein_coding")
# deg_results$results_all <- deg_results$results_all %>% filter(is.finite(log2FoldChange))


plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = gene_set_list,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_our_data_strict,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)

plot_fGSEA(deg_results = deg_results$results_all, 
           gene_set_list = sig_list_groeningen,
           title_prefix = deg_results$de_details$test, 
           save_dir = deg_dir, 
           maxSize = 3000,
           figure_extention = ".pdf"
)



vsd_counts_matrix <- assay(JunDN_vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)

annotation_col <- metadata_heatmap %>%
  dplyr::select(replicate, group) %>% 
  dplyr::arrange(group, replicate)

ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines[, match(rownames(annotation_col), colnames(ssgsea_mes_adr_cellines))]

ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        annotation_col = annotation_col,
                                                        annotation_colors = ann_colors,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = FALSE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)

pdf(file = file.path(deg_dir, "ssgsea_JunDN_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap
dev.off()


