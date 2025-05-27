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
# This is necessary to rename columns so that extraction of the data works correctly
colnames(annotationData)[c(5, 7)] <- c("ensembl_id", "gene_symbol")

OG_counts <- read.delim("~/workspace/neuroblastoma/data/RNAseq/OE_test_data/rnaseq_deseq_global_counts_raw.tsv")
OE_counts <- read.delim("~/workspace/neuroblastoma/data/RNAseq/OE_test_data/rnaseq_deseq_oe_counts_raw.tsv")
  
row.names(OG_counts) <- OG_counts$gene_id
row.names(OE_counts) <- OE_counts$gene_id

OG_counts <- OG_counts %>% select(SH_A_RNA_S130494, 
                                  SH_M_RNA_S130497)

OE_counts <- OE_counts %>% select(OE_WWTR1_DMSO_RNA_S165191, 
                                  OE_WWTR1_TDI_RNA_S165188, 
                                  OE_WWTR1_dox_RNA_S165189, 
                                  OE_WWTR1_dox_TDI_RNA_S165178,
                                  OE_mCherrry_TDI_RNA_S165184, 
                                  OE_mCherry_DMSO_RNA_S165187, 
                                  OE_mCherry_dox_RNA_S165185, 
                                  OE_mCherry_dox_TDI_RNA_S165190)

if(all(row.names(OE_counts) == row.names(OG_counts))){
  OGOE_counts <- cbind(OG_counts, OE_counts)
}
OGOE_counts
row.names(OGOE_counts) <- sub(pattern = "\\..*", 
                              replacement = "", 
                              x = row.names(OGOE_counts))

#### analysis of TDI effect WWTR_TDI/mCherry_TDI vs DMSO samples (together)
TDI_counts <- OGOE_counts %>% select(OE_WWTR1_DMSO_RNA_S165191,
                                     OE_WWTR1_TDI_RNA_S165188,
                                     OE_mCherrry_TDI_RNA_S165184,
                                     OE_mCherry_DMSO_RNA_S165187)
coldata <- data.frame(sample = c("WWTR", "WWTR", "mCherry", "mCherry"),
                      experiment = c("DMSO", "TDI", "DMSO", "TDI"))
row.names(coldata) <- colnames(TDI_counts)

TDI_dds <- DESeqDataSetFromMatrix(countData = TDI_counts,
                                  colData = coldata,
                                  design = as.formula("~ experiment"))
TDI_dds$experiment <- factor(TDI_dds$experiment, levels = c("DMSO", "TDI"))


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
  cond_interest_varPart = c("sample", "experiment"),
  color_variable = "sample",
  shape_variable = "experiment",
  ntop_genes = 1000
) +
  scale_color_manual(values = c("navy", "firebrick3")) +
  scale_shape_manual(values = c(16, 7)) +
  ggtitle("RNA-seq TDI dataset")

resultsNames(TDI_dds_filt)
deg_results <- generateResults(
  dds_object = TDI_dds_filt,
  coeff_name = "experiment_TDI_vs_DMSO",
  cond_numerator = "TDI",
  cond_denominator = "DMSO",
  cond_variable = "experiment",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff
)

# Comment - the result is non significant, the variability is too high
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
