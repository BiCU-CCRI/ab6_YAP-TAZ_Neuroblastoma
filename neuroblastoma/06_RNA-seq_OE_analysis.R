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
annotationData <- read.table(file = "~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_annotation_gene.tsv",
                             sep = "\t",
                             header = TRUE)

param_list <- list(
  abs_filt_samples=2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 0.58,
  var_expl_needed = 0.6,
  biomart_host="http://www.ensembl.org", 
  biomart_dataset="hsapiens_gene_ensembl", 
  biomart_Ens_version="Ensembl Genes 109"
)

deg_dir <- "~/workspace/neuroblastoma/results/20240229/"

#loading the new DESeq2 file
load("~/workspace/neuroblastoma/data_soren/YAP_TAZ_OverExpression/nf_out/star_salmon/deseq2_qc/deseq2.dds.RData")
new_dds_assay <- as.data.frame(counts(dds))
new_dds_assay$names_of_rows <- row.names(new_dds_assay)

wt_dds <- readRDS("~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_deseq_data_set.rds")
wt_dds_assay <- as.data.frame(counts(wt_dds))
wt_dds_assay$names_of_rows <- row.names(wt_dds_assay)

combined_counts <- merge.data.frame(new_dds_assay, 
                                    wt_dds_assay, 
                                    by = "names_of_rows") %>%
  select(-str_subset(colnames(.), "^NB*")) %>%
  rename_with(~ gsub("CLBM", "CM", .x), matches("CLBM"))
row.names(combined_counts) <- combined_counts$names_of_rows
combined_counts <- select(combined_counts, -names_of_rows)

combined_counts_col_data <- data_frame(sample_name = colnames(combined_counts))
combined_counts_col_data$cell_line <- case_when(
  str_detect(combined_counts_col_data$sample_name, "CM") ~ "CM",
  str_detect(combined_counts_col_data$sample_name, "SH") ~ "SH"
)
combined_counts_col_data$exp_type <- case_when(
  str_detect(combined_counts_col_data$sample_name, "(WWTR1)|(YAP1)") ~ "OE",
  TRUE ~ "EV"
)
combined_counts_col_data$OE_factor <- case_when(
  str_detect(combined_counts_col_data$sample_name, "WWTR1") ~ "WWTR1",
  str_detect(combined_counts_col_data$sample_name, "YAP1") ~ "YAP1",
  TRUE ~ "EV"
)
combined_counts_col_data$cell_type <- case_when(
  str_detect(combined_counts_col_data$sample_name, "_A_") ~ "ADRN",
  str_detect(combined_counts_col_data$sample_name, "_M_") ~ "MES",
  TRUE ~ "undefined"
)
combined_counts_col_data$experiment <- case_when(
  str_detect(combined_counts_col_data$sample_name, "OE") ~ "new",
  TRUE ~ "old"
)

dds <- DESeqDataSetFromMatrix(countData = combined_counts, 
                              colData = combined_counts_col_data,
                              design = ~cell_line + exp_type + cell_line:exp_type)
dds <- filterDatasets(dds, 
                      abs_filt = TRUE, 
                      abs_filt_samples = param_list$abs_filt_samples)
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::DESeq(dds)
#stabilize variance
vsd <- DESeq2::vst(dds, blind = TRUE) # blind = TRUE for QC
#generate PCA plots
pca_deg <- generatePCA_repel(transf_object = vsd, 
                             cond_interest_varPart = c("cell_type","cell_line"), 
                             color_variable = "cell_type", 
                             shape_variable = "cell_line",
                             ntop_genes = 1000) +
  ggtitle("Original dataset") 

resultsNames(dds)


# Batch correction
transf_batch_NObatch_experiment <- vsd
transf_batch_NObatch_experiment_count <- limma::removeBatchEffect(SummarizedExperiment::assay(transf_batch_NObatch_experiment), transf_batch_NObatch_experiment$experiment)
SummarizedExperiment::assay(transf_batch_NObatch_experiment) <- transf_batch_NObatch_experiment_count
pca_deg_NObatch_experiment <- generatePCA_repel(transf_object = transf_batch_NObatch_experiment, 
                                                cond_interest_varPart = c("cell_type","cell_line"), 
                                                color_variable = "cell_type", 
                                                shape_variable = "cell_line",
                                                ntop_genes = 1000) +
  ggtitle("Original dataset + YAP1/WWTR1 OE") 

pdf(file = file.path(deg_dir, "pca_wit_corrected_batch_effect.pdf"))
pca_deg_nobatch_experiment
dev.off()


#### ssGSEA analysis
# loading signatures from Groeningen paper
SIGNATURES_RAW <- READR::READ_TSV(FILE = "/HOME/RSTUDIO/WORKSPACE/NEUROBLASTOMA/MES_ADR_NCC_NORADR_GENESETS.TSV", SKIP_EMPTY_ROWS=TRUE) # DO NOT ADD NA FOR EMPTY
sig_list_groeningen <- list(
  #NCC_like = signatures_raw$`NCC-like`[!is.na(signatures_raw$`NCC-like`)],
  #Noradrenergic = signatures_raw$Noradrenergic[!is.na(signatures_raw$Noradrenergic)],
  Mesenchymal  = signatures_raw$Mesenchymal[!is.na(signatures_raw$Mesenchymal)],
  Adrenergic = signatures_raw$Adrenergic[!is.na(signatures_raw$Adrenergic)]
)
lapply(sig_list_groeningen, length)

# load the DE gene list from a previous run and assing genes to ADRN or MES identity 
# TODO ?????? check that the signatures are correct
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/20230703/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
sig_list_our_data <- list(
  Aderenergic = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "ADRN"],
  Mesenchymal = mes_adrn_gene_list$gene_symbol[mes_adrn_gene_list$Term == "MES"]
)

# prepare vsd
vsd_counts_matrix <- assay(vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_name'][match(row.names(vsd_counts_matrix), annotationData[, 'gene_id'])]

# Plot the Groeningen terms
ssgsea_mes_adr_cellines <- GSVA::gsva(vsd_counts_matrix,
                                      sig_list_groeningen,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_cellines,
                                                        scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = TRUE,
                                                        cluster_cols = TRUE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE, 
                                                        main = "SSGSEA Groeningen terms")
pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_GROEN.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap
dev.off()

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
pdf(file = file.path(deg_dir, "ssgsea_mes_adr_ncc_noradr_heatmap_OUR_RNASEQ.pdf"))
ssgsea_mes_adr_ncc_noradr_heatmap
dev.off()