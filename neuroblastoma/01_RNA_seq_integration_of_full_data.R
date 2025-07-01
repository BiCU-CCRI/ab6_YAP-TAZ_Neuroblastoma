# import::from(openxlsx, createWorkbook, addWorksheet, writeData, saveWorkbook)
# import::from(
#   .from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"),
#   "filterDatasets",
#   "generatePCA",
#   "generateEnsemblAnnotation",
#   "generateResults",
#   "meanExprsPerGroup",
#   "plotVolcano",
#   .character_only = TRUE
# ) # used for filtering
# 
# library(stringr)
# library(org.Hs.eg.db)
# library(fgsea)
# library(ggplot2)
# library(DESeq2)
# library(dplyr)
# 
# 
# # load gene sets
# gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean = TRUE)
# gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean = TRUE)
# gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
# gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean = TRUE)
# gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean = TRUE)
# gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean = TRUE)
# gs_C6_onco <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C6"), clean = TRUE)
# gs_wang_hippo <- list(wang_hippo_set = c("CCN1", "CCN2", "AMOTL2", "ANKRD1", "IGFBP3", "F3", "FJX1", "NUAK2", "LATS2", "CRIM1", "GADD45A",
#                                          "TGFB2", "PTPN14", "NT5E", "FOXF2", "AXL", "DOCK5", "ASAP1", "RBMS3", "MYOF", "ARHGEF17", "CCDC80"))
# 
# 
# 
# # Set up the parameters
# param_list <- list(
#   abs_filt_samples = 2,
#   padj_cutoff = 0.05,
#   log2FC_cutoff = 1,
#   var_expl_needed = 0.6
# )
# deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_TDI/"
# 
# annotationData <- read.table(
#   file = file.path("~/workspace/neuroblastoma/data/RNAseq/rnaseq_deseq_global_annotation_gene.tsv"),
#   sep = "\t",
#   header = TRUE
# )
# # This is necessary to rename columns so that extraction of the data works correctly
# colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")
# 
# OG_counts <- read.delim("~/workspace/neuroblastoma/data/RNAseq/OE_test_data/rnaseq_deseq_global_counts_raw.tsv")
# OE_counts <- read.delim("~/workspace/neuroblastoma/data/RNAseq/OE_test_data/rnaseq_deseq_oe_counts_raw.tsv")
#   
# row.names(OG_counts) <- OG_counts$gene_id
# row.names(OE_counts) <- OE_counts$gene_id
# 
# OG_counts <- OG_counts %>% select(SH_A_RNA_S130494, 
#                                   SH_M_RNA_S130497)
# 
# OE_counts <- OE_counts %>% select(OE_WWTR1_DMSO_RNA_S165191, 
#                                   OE_WWTR1_TDI_RNA_S165188, 
#                                   OE_WWTR1_dox_RNA_S165189, 
#                                   OE_WWTR1_dox_TDI_RNA_S165178,
#                                   OE_mCherrry_TDI_RNA_S165184, 
#                                   OE_mCherry_DMSO_RNA_S165187, 
#                                   OE_mCherry_dox_RNA_S165185, 
#                                   OE_mCherry_dox_TDI_RNA_S165190)
# 
# if(all(row.names(OE_counts) == row.names(OG_counts))){
#   OGOE_counts <- cbind(OG_counts, OE_counts)
# }
# OGOE_counts
# row.names(OGOE_counts) <- sub(pattern = "\\..*", 
#                               replacement = "", 
#                               x = row.names(OGOE_counts))
# 
# #### analysis of TDI effect WWTR_TDI/mCherry_TDI vs DMSO samples (together)
# TDI_counts <- OGOE_counts %>% select(OE_WWTR1_DMSO_RNA_S165191,
#                                      OE_WWTR1_TDI_RNA_S165188,
#                                      OE_mCherrry_TDI_RNA_S165184,
#                                      OE_mCherry_DMSO_RNA_S165187)
# coldata <- data.frame(sample = c("WWTR", "WWTR", "mCherry", "mCherry"),
#                       experiment = c("DMSO", "TDI", "DMSO", "TDI"))
# row.names(coldata) <- colnames(TDI_counts)
# 
# TDI_dds <- DESeqDataSetFromMatrix(countData = TDI_counts,
#                                   colData = coldata,
#                                   design = as.formula("~ experiment"))
# TDI_dds$experiment <- factor(TDI_dds$experiment, levels = c("DMSO", "TDI"))
# 
# 
# # make PCA with ALL samples
# # Filtering lowly expressed genes
# TDI_dds_filt <- filterDatasets(TDI_dds,
#                                       abs_filt = TRUE,
#                                       abs_filt_samples = param_list$abs_filt_samples
# )
# # Estimate size factors and running DESEQ2
# TDI_dds_filt <- DESeq2::estimateSizeFactors(TDI_dds_filt)
# TDI_dds_filt <- DESeq2::DESeq(TDI_dds_filt)
# resultsNames(TDI_dds_filt)
# # Stabilize variance
# TDI_dds_filt_vsd <- DESeq2::vst(TDI_dds_filt, blind = TRUE) # blind = TRUE for QC
# 
# # generate PCA plots
# tmp_pca_deg <- generatePCA(
#   transf_object = TDI_dds_filt_vsd,
#   cond_interest_varPart = c("sample", "experiment"),
#   color_variable = "sample",
#   shape_variable = "experiment",
#   ntop_genes = 1000
# ) +
#   scale_color_manual(values = c("navy", "firebrick3")) +
#   scale_shape_manual(values = c(16, 7)) +
#   ggtitle("RNA-seq TDI dataset")
# 
# resultsNames(TDI_dds_filt)
# deg_results <- generateResults(
#   dds_object = TDI_dds_filt,
#   coeff_name = "experiment_TDI_vs_DMSO",
#   cond_numerator = "TDI",
#   cond_denominator = "DMSO",
#   cond_variable = "experiment",
#   ensemblAnnot = annotationData,
#   log2FC_cutoff = param_list$log2FC_cutoff
# )
# 
# # Comment - the result is non significant, the variability is too high
# gene_set_list <- list(
#   GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
#   REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
#   C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
#   Hippo_Wang = gs_wang_hippo$wang_hippo_set
# )
# 
# de_genes_ranked <- deg_results$results_all
# de_genes_ranked <- de_genes_ranked %>%
#   group_by(gene_symbol) %>%
#   slice_max(baseMean, n = 1, with_ties = FALSE) %>%
#   ungroup()
# de_genes_ranked <- de_genes_ranked %>% arrange(desc(log2FoldChange))
# de_genes_ranked <- setNames(c(de_genes_ranked$log2FoldChange), c(de_genes_ranked$gene_symbol))
# 
# 
# fgseaRes <- fgsea(
#   pathways = gene_set_list,
#   stats = de_genes_ranked,
#   minSize = 15,
#   maxSize = 500
# )
# 
# for(gene_set_to_plot_name in names(gene_set_list)){
#   # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
#   gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
#   
#   title_values <- fgseaRes %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
#   p <- plotEnrichment(pathway = gene_set_to_plot, stats = de_genes_ranked) + 
#     labs(title = paste(title_values$pathway, "\n",
#                        "NES =", title_values$NES, "\n",
#                        "padj = ", title_values$padj))
#   plot(p)
#   
#   ggsave(
#     filename = paste0(deg_dir, "fgsea_",gene_set_to_plot_name, ".pdf"),
#     plot = p,
#     width = 20, height = 20, units = "cm"
#   )
# }
# 
# 
# 

###########################################################################################################################
## new analysis whith not united replicates 
import::from(openxlsx, createWorkbook, addWorksheet, writeData, saveWorkbook)
import::from(
  .from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"),
  "filterDatasets",
  "generatePCA",
  "generateEnsemblAnnotation",
  "generateResults",
  "meanExprsPerGroup",
  "plotVolcano",
  .character_only = TRUE
) # used for filtering

library(stringr)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(DESeq2)
library(dplyr)


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



# Set up the parameters
param_list <- list(
  abs_filt_samples = 2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 1,
  var_expl_needed = 0.6
)
deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_TDI/"

annotationData <- read.table(
  file = file.path("~/workspace/neuroblastoma/data/RNAseq/rnaseq_deseq_global_annotation_gene.tsv"),
  sep = "\t",
  header = TRUE
)
colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")

load("~/workspace/neuroblastoma/data/RNAseq/OE_WWTR1_mCherry_JunDN_experiment.dds.RData")
samples_to_select_full <- c("OE_WWTR1_DMSO_RNA_S165191_R1", 
                            "OE_WWTR1_TDI_RNA_S165188_R1", 
                            "OE_WWTR1_dox_RNA_S165189_R1", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R1",
                            "OE_mCherrry_TDI_RNA_S165184_R1", 
                            "OE_mCherry_DMSO_RNA_S165187_R1", 
                            "OE_mCherry_dox_RNA_S165185_R1", 
                            "OE_mCherry_dox_TDI_RNA_S165190_R1",
                            "OE_WWTR1_DMSO_RNA_S165191_R2", 
                            "OE_WWTR1_TDI_RNA_S165188_R2", 
                            "OE_WWTR1_dox_RNA_S165189_R2", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R2",
                            "OE_mCherrry_TDI_RNA_S165184_R2", 
                            "OE_mCherry_DMSO_RNA_S165187_R2", 
                            "OE_mCherry_dox_RNA_S165185_R2", 
                            "OE_mCherry_dox_TDI_RNA_S165190_R2")

samples_to_select_TDI_vs_DMSO <- c("OE_WWTR1_DMSO_RNA_S165191_R1", 
                                   "OE_WWTR1_TDI_RNA_S165188_R1", 
                                   "OE_mCherrry_TDI_RNA_S165184_R1", 
                                   "OE_mCherry_DMSO_RNA_S165187_R1", 
                                   "OE_WWTR1_DMSO_RNA_S165191_R2", 
                                   "OE_WWTR1_TDI_RNA_S165188_R2", 
                                   "OE_mCherrry_TDI_RNA_S165184_R2", 
                                   "OE_mCherry_DMSO_RNA_S165187_R2")

TDI_dds <- dds[,colnames(dds) %in% samples_to_select_TDI_vs_DMSO]
#TDI_dds <- dds[,colnames(dds) %in% samples_to_select_full]

sample_names <- colData(TDI_dds)$sample

parsed <- str_match(
  sample_names,
  "^OE_([^_]+)_((?:[^_]+_)*[^_]+)_RNA_.*_(R\\d+)$"
)
parsed[,2][parsed[, 2] == "mCherrry"] <- "mCherry"

# Create a tidy data frame
parsed_df <- data.frame(
  genotype  = parsed[, 2],
  treatment = parsed[, 3],
  replicate = parsed[, 4],
  stringsAsFactors = TRUE
)

colData(TDI_dds) <- cbind(colData(TDI_dds), parsed_df)
#TDI_dds$treatment <- factor(TDI_dds$treatment, levels = c("DMSO", "TDI", "dox", "dox_TDI"))
TDI_dds$treatment <- factor(TDI_dds$treatment, levels = c("DMSO", "TDI"))

design(TDI_dds) <- as.formula("~ genotype + treatment")



# make PCA with ALL samples
# Filtering lowly expressed genes
TDI_dds_filt <- filterDatasets(TDI_dds,
                               abs_filt = TRUE,
                               abs_filt_samples = param_list$abs_filt_samples
)

# Estimate size factors and running DESEQ2
TDI_dds_filt <- DESeq2::estimateSizeFactors(TDI_dds_filt)
TDI_dds_filt <- DESeq2::DESeq(TDI_dds_filt)
resultsNames(TDI_dds_filt)
# Stabilize variance
TDI_dds_filt_vsd <- DESeq2::vst(TDI_dds_filt, blind = TRUE) # blind = TRUE for QC

# generate PCA plots
tmp_pca_deg <- generatePCA(
  transf_object = TDI_dds_filt_vsd,
  cond_interest_varPart = c("genotype", "treatment"),
  color_variable = "genotype",
  shape_variable = "treatment",
  ntop_genes = 1000
) +
  scale_color_manual(values = c("navy", "firebrick3")) +
  scale_shape_manual(values = c(16, 7)) +
  ggtitle("RNA-seq TDI dataset")

resultsNames(TDI_dds_filt)

deg_results <- generateResults(
  dds_object = TDI_dds_filt,
  coeff_name = "treatment_TDI_vs_DMSO",
  cond_numerator = "TDI",
  cond_denominator = "DMSO",
  cond_variable = "treatment",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

deg_results_wwtr_mcherry <- deg_results

####
gene_set_list <- list(
  GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
  REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
  C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
  Hippo_Wang = gs_wang_hippo$wang_hippo_set
)

de_genes_ranked <- deg_results$results_all
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
  maxSize = 500
)

for(gene_set_to_plot_name in names(gene_set_list)){
  # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
  gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
  
  title_values <- fgseaRes %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
  p <- plotEnrichment(pathway = gene_set_to_plot, stats = de_genes_ranked) + 
    labs(title = paste(title_values$pathway, "\n",
                       "NES =", title_values$NES, "\n",
                       "padj = ", title_values$padj))
  plot(p)
  
  ggsave(
    filename = paste0(deg_dir, "fgsea_",gene_set_to_plot_name, ".pdf"),
    plot = p,
    width = 20, height = 20, units = "cm"
  )
}


#####
# GSVA using our terms
# load the DE gene list from a previous run and assing genes to ADRN or MES identity 
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)

# prepare vsd
vsd_counts_matrix <- assay(TDI_dds_filt_vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

# Plot terms based on terms defined from RNA-seq profiles obtained from our RNA seq
ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = TRUE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)
#pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap




########################################### 
# separate mCherry and WWTR
deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_TDI/"

annotationData <- read.table(
  file = file.path("~/workspace/neuroblastoma/data/RNAseq/rnaseq_deseq_global_annotation_gene.tsv"),
  sep = "\t",
  header = TRUE
)
colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")

load("~/workspace/neuroblastoma/data/RNAseq/OE_WWTR1_mCherry_JunDN_experiment.dds.RData")

samples_to_select_TDI_vs_DMSO <- c("OE_WWTR1_DMSO_RNA_S165191_R1", 
                                   "OE_WWTR1_TDI_RNA_S165188_R1", 
                                   "OE_WWTR1_DMSO_RNA_S165191_R2", 
                                   "OE_WWTR1_TDI_RNA_S165188_R2"
                                   )

TDI_dds <- dds[,colnames(dds) %in% samples_to_select_TDI_vs_DMSO]
sample_names <- colData(TDI_dds)$sample

parsed <- str_match(
  sample_names,
  "^OE_([^_]+)_((?:[^_]+_)*[^_]+)_RNA_.*_(R\\d+)$"
)

# Create a tidy data frame
parsed_df <- data.frame(
  genotype  = parsed[, 2],
  treatment = parsed[, 3],
  replicate = parsed[, 4],
  stringsAsFactors = TRUE
)

colData(TDI_dds) <- cbind(colData(TDI_dds), parsed_df)
#TDI_dds$treatment <- factor(TDI_dds$treatment, levels = c("DMSO", "TDI", "dox", "dox_TDI"))
TDI_dds$treatment <- factor(TDI_dds$treatment, levels = c("DMSO", "TDI"))

design(TDI_dds) <- as.formula("~ treatment")

# make PCA with ALL samples
# Filtering lowly expressed genes
TDI_dds_filt <- filterDatasets(TDI_dds,
                               abs_filt = TRUE,
                               abs_filt_samples = param_list$abs_filt_samples
)

# Estimate size factors and running DESEQ2
TDI_dds_filt <- DESeq2::estimateSizeFactors(TDI_dds_filt)
TDI_dds_filt <- DESeq2::DESeq(TDI_dds_filt)
resultsNames(TDI_dds_filt)
# Stabilize variance
TDI_dds_filt_vsd <- DESeq2::vst(TDI_dds_filt, blind = TRUE) # blind = TRUE for QC

# generate PCA plots
tmp_pca_deg <- generatePCA(
  transf_object = TDI_dds_filt_vsd,
  cond_interest_varPart = c("genotype", "treatment"),
  color_variable = "genotype",
  shape_variable = "treatment",
  ntop_genes = 1000
) +
  scale_color_manual(values = c("navy", "firebrick3")) +
  scale_shape_manual(values = c(16, 7)) +
  ggtitle("RNA-seq TDI dataset")

resultsNames(TDI_dds_filt)

deg_results <- generateResults(
  dds_object = TDI_dds_filt,
  coeff_name = "treatment_TDI_vs_DMSO",
  cond_numerator = "TDI",
  cond_denominator = "DMSO",
  cond_variable = "treatment",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

####
gene_set_list <- list(
  GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
  REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
  C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
  Hippo_Wang = gs_wang_hippo$wang_hippo_set
)

de_genes_ranked <- deg_results$results_all
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
  maxSize = 500
)

for(gene_set_to_plot_name in names(gene_set_list)){
  # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
  gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
  
  title_values <- fgseaRes %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
  p <- plotEnrichment(pathway = gene_set_to_plot, stats = de_genes_ranked) + 
    labs(title = paste(title_values$pathway, "\n",
                       "NES =", title_values$NES, "\n",
                       "padj = ", title_values$padj))
  plot(p)
  
  ggsave(
    filename = paste0(deg_dir, "WWTR_only_fgsea_",gene_set_to_plot_name, ".pdf"),
    plot = p,
    width = 20, height = 20, units = "cm"
  )
}


#####
# GSVA using our terms
# load the DE gene list from a previous run and assing genes to ADRN or MES identity 
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)

# prepare vsd
vsd_counts_matrix <- assay(TDI_dds_filt_vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

# Plot terms based on terms defined from RNA-seq profiles obtained from our RNA seq
ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = TRUE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)
#pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap





#### TEST compare the WWTR only and WWTR mCherry experiment
deg_results$results_signif
write.table(
  x = deg_results$results_signif,
  file = "~/workspace/neuroblastoma/temp_results/deg_results_wwtr.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)


write.table(
  x = deg_results_wwtr_mcherry$results_signif,
  file = "~/workspace/neuroblastoma/temp_results/deg_results_wwtr_mcherry.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)

intersect(deg_results$results_signif$gene_symbol, 
          deg_results_wwtr_mcherry$results_signif$gene_symbol)











###########################################################################################################################
## new analysis whith not united replicates TDI, dox, dox_TDI, DMSO
import::from(openxlsx, createWorkbook, addWorksheet, writeData, saveWorkbook)
import::from(
  .from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"),
  "filterDatasets",
  "generatePCA",
  "generateEnsemblAnnotation",
  "generateResults",
  "meanExprsPerGroup",
  "plotVolcano",
  .character_only = TRUE
) # used for filtering

library(stringr)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(DESeq2)
library(dplyr)


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



# Set up the parameters
param_list <- list(
  abs_filt_samples = 2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 1,
  var_expl_needed = 0.6
)
deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq_TDI/"

annotationData <- read.table(
  file = file.path("~/workspace/neuroblastoma/data/RNAseq/rnaseq_deseq_global_annotation_gene.tsv"),
  sep = "\t",
  header = TRUE
)
colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")

load("~/workspace/neuroblastoma/data/RNAseq/OE_WWTR1_mCherry_JunDN_experiment.dds.RData")
samples_to_select_full <- c("OE_WWTR1_DMSO_RNA_S165191_R1", 
                            "OE_WWTR1_TDI_RNA_S165188_R1", 
                            "OE_WWTR1_dox_RNA_S165189_R1", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R1",
                            "OE_WWTR1_DMSO_RNA_S165191_R2", 
                            "OE_WWTR1_TDI_RNA_S165188_R2", 
                            "OE_WWTR1_dox_RNA_S165189_R2", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R2"
                            )

TDI_dds <- dds[,colnames(dds) %in% samples_to_select_full]

sample_names <- colData(TDI_dds)$sample

parsed <- str_match(
  sample_names,
  "^OE_([^_]+)_((?:[^_]+_)*[^_]+)_RNA_.*_(R\\d+)$"
)

# Create a tidy data frame
parsed_df <- data.frame(
  genotype  = parsed[, 2],
  treatment = parsed[, 3],
  replicate = parsed[, 4],
  stringsAsFactors = TRUE
)

colData(TDI_dds) <- cbind(colData(TDI_dds), parsed_df)
TDI_dds$treatment <- factor(TDI_dds$treatment, levels = c("DMSO", "TDI", "dox", "dox_TDI"))

design(TDI_dds) <- as.formula("~ treatment")

# make PCA with ALL samples
# Filtering lowly expressed genes
TDI_dds_filt <- filterDatasets(TDI_dds,
                               abs_filt = TRUE,
                               abs_filt_samples = param_list$abs_filt_samples
)

# Estimate size factors and running DESEQ2
TDI_dds_filt <- DESeq2::estimateSizeFactors(TDI_dds_filt)
TDI_dds_filt <- DESeq2::DESeq(TDI_dds_filt)
resultsNames(TDI_dds_filt)

deg_results_TDI_vs_DMSO <- generateResults(
  dds_object = TDI_dds_filt,
  coeff_name = "treatment_TDI_vs_DMSO",
  cond_numerator = "TDI",
  cond_denominator = "DMSO",
  cond_variable = "treatment",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

# Write the result to the xlsx file
XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")
openxlsx::writeData(XLSX_OUT, x = deg_results_TDI_vs_DMSO$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = deg_results_TDI_vs_DMSO$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = deg_results_TDI_vs_DMSO$results_all, sheet = "results_all")
openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir, "treatment_TDI_vs_DMSO.xlsx"), overwrite = T)

deg_results_doxTDI_vs_dox <- generateResults(
  dds_object = TDI_dds_filt,
  coef_contrast=c("treatment", "dox_TDI", "dox"),
  cond_numerator = "dox_TDI",
  cond_denominator = "dox",
  cond_variable = "treatment",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

# Write the result to the xlsx file
XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")

openxlsx::writeData(XLSX_OUT, x = deg_results_doxTDI_vs_dox$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = deg_results_doxTDI_vs_dox$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = deg_results_doxTDI_vs_dox$results_all, sheet = "results_all")

openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir, "treatment_doxTDI_vs_dox.xlsx.xlsx"), overwrite = T)

# Stabilize variance
TDI_dds_filt_vsd <- DESeq2::vst(TDI_dds_filt, blind = TRUE) # blind = TRUE for QC

#remove R1 replicates
# generate PCA plots
tmp_pca_deg <- generatePCA(
  transf_object = TDI_dds_filt_vsd[, c(1,3,5,7)],
  cond_interest_varPart = c("genotype", "treatment"),
  color_variable = "genotype",
  shape_variable = "treatment",
  ntop_genes = 1000
) + ggtitle("PCA plot for dox doxTDI DMSO TDI WWTR1 only")
ggsave(
  filename = paste0(deg_dir, "PCA_WWTR1_only_dox_doxTDI_DMSO.pdf"),
  plot = tmp_pca_deg,
  width = 20, height = 20, units = "cm"
)

#####################################################
gene_set_list <- list(
  GOBP_Hippo_Signaling = gs_C5_GOBP[["genesets"]][["Hippo Signaling"]],
  REACTOME_Signaling_By_Hippo = gs_C2_reactome[["genesets"]][["Signaling By Hippo"]],
  C6_onko_Cordenonsi_Yap_Conserved_Signature = gs_C6_onco[["genesets"]][["Cordenonsi Yap Conserved Signature"]],
  Hippo_Wang = gs_wang_hippo$wang_hippo_set
)

# TDI_vs_DMSO
deg_results_TDI_vs_DMSO_ranked <- deg_results_TDI_vs_DMSO$results_all
deg_results_TDI_vs_DMSO_ranked <- deg_results_TDI_vs_DMSO_ranked %>%
  group_by(gene_symbol) %>%
  slice_max(baseMean, n = 1, with_ties = FALSE) %>%
  ungroup()
deg_results_TDI_vs_DMSO_ranked <- deg_results_TDI_vs_DMSO_ranked %>% arrange(desc(log2FoldChange))
deg_results_TDI_vs_DMSO_ranked <- setNames(c(deg_results_TDI_vs_DMSO_ranked$log2FoldChange), c(deg_results_TDI_vs_DMSO_ranked$gene_symbol))


fgseaRes_TDI_vs_DMSO <- fgsea(
  pathways = gene_set_list,
  stats = deg_results_TDI_vs_DMSO_ranked,
  minSize = 15,
  maxSize = 500
)

for(gene_set_to_plot_name in names(gene_set_list)){
  # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
  gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
  
  title_values <- fgseaRes_TDI_vs_DMSO %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
  p <- plotEnrichment(pathway = gene_set_to_plot, stats = deg_results_TDI_vs_DMSO_ranked) + 
    labs(title = paste(title_values$pathway, "\n",
                       "NES =", title_values$NES, "\n",
                       "padj = ", title_values$padj))
  plot(p)
  
  ggsave(
    filename = paste0(deg_dir, "TDI_vs_DMSO_fgsea_", gene_set_to_plot_name, ".pdf"),
    plot = p,
    width = 20, height = 20, units = "cm"
  )
}

############## doxTDI_vs_dox
deg_results_doxTDI_vs_dox_ranked <- deg_results_doxTDI_vs_dox$results_all
deg_results_doxTDI_vs_dox_ranked <- deg_results_doxTDI_vs_dox_ranked %>%
  group_by(gene_symbol) %>%
  slice_max(baseMean, n = 1, with_ties = FALSE) %>%
  ungroup()
deg_results_doxTDI_vs_dox_ranked <- deg_results_doxTDI_vs_dox_ranked %>% arrange(desc(log2FoldChange))
deg_results_doxTDI_vs_dox_ranked <- setNames(c(deg_results_doxTDI_vs_dox_ranked$log2FoldChange), c(deg_results_doxTDI_vs_dox_ranked$gene_symbol))


fgseaRes_doxTDI_vs_dox <- fgsea(
  pathways = gene_set_list,
  stats = deg_results_doxTDI_vs_dox_ranked,
  minSize = 15,
  maxSize = 500
)

for(gene_set_to_plot_name in names(gene_set_list)){
  # gene_set_to_plot_name <- "Cordenonsi_Yap_Conserved_Signature"
  gene_set_to_plot <- gene_set_list[[gene_set_to_plot_name]]
  
  title_values <- fgseaRes_doxTDI_vs_dox %>% filter(pathway == gene_set_to_plot_name) %>% select(pathway, NES, padj)
  p <- plotEnrichment(pathway = gene_set_to_plot, stats = deg_results_doxTDI_vs_dox_ranked) + 
    labs(title = paste(title_values$pathway, "\n",
                       "NES =", title_values$NES, "\n",
                       "padj = ", title_values$padj))
  plot(p)
  
  ggsave(
    filename = paste0(deg_dir, "doxTDI_vs_dox_fgsea_", gene_set_to_plot_name, ".pdf"),
    plot = p,
    width = 20, height = 20, units = "cm"
  )
}

################################################################################
# Add old ADR and MES samples for PCA visualisation
load("~/workspace/neuroblastoma/data/RNAseq/OE_WWTR1_mCherry_JunDN_experiment.dds.RData")
samples_to_select_full <- c("OE_WWTR1_DMSO_RNA_S165191_R1", 
                            "OE_WWTR1_TDI_RNA_S165188_R1", 
                            "OE_WWTR1_dox_RNA_S165189_R1", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R1",
                            "OE_mCherrry_TDI_RNA_S165184_R1", 
                            "OE_mCherry_DMSO_RNA_S165187_R1", 
                            "OE_mCherry_dox_RNA_S165185_R1", 
                            "OE_mCherry_dox_TDI_RNA_S165190_R1",
                            "OE_WWTR1_DMSO_RNA_S165191_R2", 
                            "OE_WWTR1_TDI_RNA_S165188_R2", 
                            "OE_WWTR1_dox_RNA_S165189_R2", 
                            "OE_WWTR1_dox_TDI_RNA_S165178_R2",
                            "OE_mCherrry_TDI_RNA_S165184_R2", 
                            "OE_mCherry_DMSO_RNA_S165187_R2", 
                            "OE_mCherry_dox_RNA_S165185_R2", 
                            "OE_mCherry_dox_TDI_RNA_S165190_R2")
TDI_dds <- dds[ ,colnames(dds) %in% samples_to_select_full]

samples_to_select_OG <- c("SH_A_RNA_S130494","SH_M_RNA_S130497")
rnaseq_OG <- readRDS("~/workspace/neuroblastoma/data/RNAseq/rnaseq_deseq_global_deseq_data_set.rds")
rnaseq_OG <- rnaseq_OG[ ,colnames(rnaseq_OG) %in% samples_to_select_OG]

# Assuming dds1 and dds2 are your DESeq2 objects
TDI_dds_counts <- as.data.frame(TDI_dds@assays@data@listData$counts)
rnaseq_OG_counts <- as.data.frame(rnaseq_OG@assays@data@listData$counts)

TDI_dds_counts$gene_id <- rownames(TDI_dds)
rnaseq_OG_counts$gene_id <- rownames(rnaseq_OG)

TDI_plus_OG_counts <- merge.data.frame(TDI_dds_counts, rnaseq_OG_counts, by = "gene_id")
row.names(TDI_plus_OG_counts) <- TDI_plus_OG_counts$gene_id
TDI_plus_OG_counts <- select(TDI_plus_OG_counts, !gene_id)

parsed <- str_match(
  colnames(TDI_plus_OG_counts),
  "^OE_([^_]+)_((?:[^_]+_)*[^_]+)_RNA_.*_(R\\d+)$"
)
parsed[,2][parsed[, 2] == "mCherrry"] <- "mCherry"

#manual additon of two control samples
parsed[c(17,18), 1] <- colnames(TDI_plus_OG_counts)[c(17,18)]
parsed[c(17,18), 2] <- c("outside_control_ADR", "outside_control_MES")
parsed[c(17,18), 3] <- c("ADR", "MES")
parsed[c(17,18), 4] <- "R1"

# Create a tidy data frame
parsed_df <- data.frame(
  genotype  = as.factor(parsed[, 2]),
  treatment = as.factor(parsed[, 3]),
  replicate = as.factor(parsed[, 4]),
  stringsAsFactors = TRUE
)
parsed_df$is_outer_control <- as.factor(c(rep("not_control", 16), rep("control", 2)))

TDI_plus_OG <- DESeqDataSetFromMatrix(countData = TDI_plus_OG_counts,
                                  colData = parsed_df,
                                  design = as.formula("~ 1"))

# design(TDI_plus_OG) <- as.formula("~1")

# make PCA with ALL samples
# Filtering lowly expressed genes


TDI_plus_OG_filt <- filterDatasets(TDI_plus_OG,
                               abs_filt = TRUE,
                               abs_filt_samples = param_list$abs_filt_samples
)

# ADD OR REMOVE mCHERRY samples
TDI_plus_OG_filt <- TDI_plus_OG_filt[, colData(TDI_plus_OG_filt)$genotype != "mCherry"]

# Estimate size factors and running DESEQ2
TDI_plus_OG_filt <- DESeq2::estimateSizeFactors(TDI_plus_OG_filt)
TDI_plus_OG_filt <- DESeq2::DESeq(TDI_plus_OG_filt)
resultsNames(TDI_plus_OG_filt)
# Stabilize variance
TDI_plus_OG_filt_vsd <- DESeq2::vst(TDI_plus_OG_filt, blind = TRUE) # blind = TRUE for QC

# generate PCA plots
tmp_pca_deg <- generatePCA(
  transf_object = TDI_plus_OG_filt_vsd,
  cond_interest_varPart = c("genotype", "treatment"),
  color_variable = "genotype",
  shape_variable = "treatment",
  ntop_genes = 1000
) +
ggtitle("dox tdi dox/tdi dmso experiment + OG samples. NOT batch corrected") 
tmp_pca_deg

# Batch correction
transf_batch_NObatch_experiment <- TDI_plus_OG_filt_vsd
transf_batch_NObatch_experiment_count <- limma::removeBatchEffect(SummarizedExperiment::assay(transf_batch_NObatch_experiment),
                                                                  transf_batch_NObatch_experiment$is_outer_control)
SummarizedExperiment::assay(transf_batch_NObatch_experiment) <- transf_batch_NObatch_experiment_count
pca_deg_NObatch_experiment <- generatePCA(transf_object = transf_batch_NObatch_experiment[,c(1,3,5,7,9,10)], 
                                                cond_interest_varPart = c("genotype","treatment"), 
                                                color_variable = "genotype", 
                                                shape_variable = "treatment",
                                                ntop_genes = 1000) +
  ggtitle("dox tdi dox/tdi dmso experiment + OG samples. batch corrected") 
pca_deg_NObatch_experiment

ggsave(
  filename = paste0(deg_dir, "PCA_WWTR1_only_dox_doxTDI_DMSO_plus_OGs.pdf"),
  plot = pca_deg_NObatch_experiment,
  width = 20, height = 20, units = "cm"
)

########## GSVA
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)

# prepare vsd
vsd_counts_matrix <- assay(TDI_plus_OG_filt_vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

# Plot terms based on terms defined from RNA-seq profiles obtained from our RNA seq
ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_our_data,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                      #  scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = FALSE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)
#pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap