###
# Title: RNA-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
#

# Set up the environment ####
# Importing only key functions that are actually used - not to polute namespace!
import::from(openxlsx, createWorkbook, addWorksheet, writeData, saveWorkbook)
import::from(
  .from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"),
  "filterDatasets",
  "generatePCA",
  "generateEnsemblAnnotation",
  "generateResults",
  "meanExprsPerGroup",
  "plotVolcano",
  .character_only = TRUE) # used for filtering

library(xcore)
#library(ExperimentHub)
#library(xcoredata)
library(stringr)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(DESeq2)

# Load 3rd party packages
# Note - might take several attempts to actually load it
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org")
# Loading MsigDB geneset collections
gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean = TRUE)
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean = TRUE)
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean = TRUE)
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean = TRUE)
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean = TRUE)

# Set up the parameters
param_list <- list(
  abs_filt_samples=2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 1,
  var_expl_needed = 0.6
  #biomart_host="http://www.ensembl.org",
  #biomart_dataset="hsapiens_gene_ensembl",
  #biomart_Ens_version="Ensembl Genes 109"
  )
deg_dir <- "~/workspace/neuroblastoma/results/RNA-seq/"



# Analysis of the RNA-seq data ####
# Loading the DESeq2 file and annotation file that was used during the pre-processing of the data
path_folder_rna_seq_data <- "~/workspace/neuroblastoma/data/" 
rnaseq_dds <- readRDS(file.path(path_folder_rna_seq_data, "rnaseq_deseq_global_deseq_data_set.rds"))
design(rnaseq_dds) <- as.formula("~ cell_line + cell_type")
annotationData <- read.table(file = file.path(path_folder_rna_seq_data, "rnaseq_deseq_global_annotation_gene.tsv"),
                             sep = "\t",
                             header = TRUE)
# This is necessary to rename columns so that extraction of the data works correctly
colnames(annotationData)[c(5,7)] <- c("ensembl_id", "gene_symbol")


# make PCA with ALL samples
# Filtering lowly expressed genes
tmp_rnaseq_dds_filt <- filterDatasets(rnaseq_dds, 
                                      abs_filt = TRUE, 
                                      abs_filt_samples = param_list$abs_filt_samples)
# Estimate size factors and running DESEQ2
tmp_rnaseq_dds_filt <- DESeq2::estimateSizeFactors(tmp_rnaseq_dds_filt)
tmp_rnaseq_dds_filt <- DESeq2::DESeq(tmp_rnaseq_dds_filt)
# Stabilize variance
tmp_vsd <- DESeq2::vst(tmp_rnaseq_dds_filt, blind = TRUE) # blind = TRUE for QC

#generate PCA plots
tmp_pca_deg <- generatePCA(transf_object = tmp_vsd, 
                           cond_interest_varPart = c("cell_type","cell_line"), 
                           color_variable = "cell_type", 
                           shape_variable = "cell_line",
                           ntop_genes = 1000) + 
  scale_color_manual(values = c("navy", "firebrick3")) +
  scale_shape_manual(values = c(16,7,17,15,3))+
  ggtitle("RNA-seq Original dataset") 
ggsave(filename = paste0(deg_dir, "PCA_full_dataset.png"), 
       plot = tmp_pca_deg,
       width = 20, height = 20, units = "cm")

# Removing NB6 sample - this is intermixed sample
rnaseq_dds <- rnaseq_dds[,rnaseq_dds$cell_line != "STA_NB_6"]
rnaseq_dds$cell_line  <- droplevels(rnaseq_dds$cell_line)
# Filtering lowly expressed genes
rnaseq_dds_filt <- filterDatasets(rnaseq_dds, 
                                  abs_filt = TRUE, 
                                  abs_filt_samples = param_list$abs_filt_samples)
# Estimate size factors and running DESEQ2
rnaseq_dds_filt <- DESeq2::estimateSizeFactors(rnaseq_dds_filt)
rnaseq_dds_filt <- DESeq2::DESeq(rnaseq_dds_filt)
# Stabilize variance
vsd <- DESeq2::vst(rnaseq_dds_filt, blind = TRUE) # blind = TRUE for QC

#generate PCA plots
pca_deg <- generatePCA(
  transf_object = vsd, 
  cond_interest_varPart = c("cell_type","cell_line"), 
  color_variable = "cell_type", 
  shape_variable = "cell_line",
  ntop_genes = 1000) +
  scale_color_manual(values = c("navy", "firebrick3")) +
  ggtitle("Original dataset")
ggsave(filename = file.path(deg_dir, "PCA_MES_vs_ADRN_DEG.png"), 
       plot = pca_deg,
       width = 20, height = 20, units = "cm")

#   
# ensemblAnnot <- generateEnsemblAnnotation(ensembl_ids = rownames(rnaseq_dds_filt),
#                                           host=param_list$biomart_host,
#                                           version=param_list$biomart_Ens_version,
#                                           dataset=param_list$biomart_dataset)

# Check the names of the dds object and generate results
resultsNames(rnaseq_dds_filt)
deg_results <- generateResults(
  dds_object = rnaseq_dds_filt, 
  coeff_name = "cell_type_MES_vs_ADR",
  cond_numerator = "MES", 
  cond_denominator = "ADR",
  cond_variable="cell_type",
  ensemblAnnot = annotationData,
  log2FC_cutoff = param_list$log2FC_cutoff)

#export MES and ADR genes and create custom terms from them
ADR_genes <- deg_results$results_signif %>% filter(log2FoldChange < 0) %>% pull(gene_symbol)
MES_genes <- deg_results$results_signif %>% filter(log2FoldChange > 0) %>% pull(gene_symbol)
hs <- org.Hs.eg.db
ADR_genes <- AnnotationDbi::select(
  hs,
  keys = ADR_genes,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL")
ADR_genes <- ADR_genes[!is.na(ADR_genes$ENTREZID),]
ADR_genes$gs_id <- "ADRN"

MES_genes <- AnnotationDbi::select(hs,
                                   keys = MES_genes,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
MES_genes <- MES_genes[!is.na(MES_genes$ENTREZID),]
MES_genes$gs_id <- "MES"

ADRN_MES_terms <- data.frame(
  rbind(
    MES_genes %>% dplyr::select(gs_id, ENTREZID),
    ADR_genes %>% dplyr::select(gs_id, ENTREZID)
  )
)
colnames(ADRN_MES_terms) <- c("gs_id",	"gene_id")
write.table(ADRN_MES_terms, 
            file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_from_RNA_seq.tsv", 
            quote = FALSE,
            sep = "\t",
            row.names = FALSE, 
            col.names = TRUE)

# Save results to an xlsx sheet
XLSX_OUT <- createWorkbook()
addWorksheet(XLSX_OUT, "results_signif")
addWorksheet(XLSX_OUT, "de_details")
addWorksheet(XLSX_OUT, "results_all")

writeData(XLSX_OUT, x = deg_results$results_signif, sheet = "results_signif")
writeData(XLSX_OUT, x = deg_results$de_details, sheet = "de_details")
writeData(XLSX_OUT, x = deg_results$results_all, sheet = "results_all")

saveWorkbook(XLSX_OUT, file.path(deg_dir, "cell_type_MES_vs_ADR.xlsx"))

# Produce GSEA plots
C5_GOBP <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = gs_C5_GOBP, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))
C5_GOCC <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = gs_C5_GOCC, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))
C5_GOMF <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = gs_C5_GOMF, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))
C2_kegg <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = gs_C2_kegg, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))
C2_reactome <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = gs_C2_reactome, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))

C5_GOBP_plot <- hypeR::hyp_dots(C5_GOBP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: MES vs ADR") +theme_bw()
C5_GOCC_plot <- hypeR::hyp_dots(C5_GOCC, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOCC: MES vs ADR") +theme_bw()
C5_GOMF_plot <- hypeR::hyp_dots(C5_GOMF, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOMF: MES vs ADR") +theme_bw()
C2_kegg_plot <- hypeR::hyp_dots(C2_kegg, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="KEGG: MES vs ADR") +theme_bw()
C2_rctm_plot <- hypeR::hyp_dots(C2_reactome, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="REACTOME: MES vs ADR") +theme_bw()

#optional saving to excel tables
#hypeR::hyp_to_excel(C5_GOBP, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOBP.xlsx"))
#hypeR::hyp_to_excel(C5_GOCC, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOCC.xlsx"))
#hypeR::hyp_to_excel(C5_GOMF, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOMF.xlsx"))

pdf(file = file.path(deg_dir, "MES_vs_ADR_GSEA.pdf"), width = 10, height = 10)
plot(C5_GOBP_plot)
plot(C5_GOCC_plot)
plot(C5_GOMF_plot)
plot(C2_kegg_plot)
plot(C2_rctm_plot)
dev.off()

# Additionally - save as individual plots
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_GSEA_C5_GOBP.png"), 
       plot = C5_GOBP_plot,
       width = 20, height = 20, units = "cm")
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_GSEA_C5_GOCC.png"), 
       plot = C5_GOCC_plot,
       width = 20, height = 20, units = "cm")
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_GSEA_C5_GOMF.png"), 
       plot = C5_GOMF_plot,
       width = 20, height = 20, units = "cm")
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_GSEA_C2_KEGG.png"), 
       plot = C2_kegg_plot,
       width = 20, height = 20, units = "cm")
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_GSEA_C2_REACTOME.png"), 
       plot = C2_rctm_plot,
       width = 20, height = 20, units = "cm")

# Run Fast GSEA for GOBP, GOCC, GOMF, KEGG, reactome terms

# Rank genes based on L2FC
de_genes_ranked <- deg_results$results_all
de_genes_ranked <- de_genes_ranked %>% arrange(desc(log2FoldChange))
de_genes_ranked <- setNames(c(de_genes_ranked$log2FoldChange), c(de_genes_ranked$gene_symbol))


gene_set_list <- list(
  gs_hallmark = gs_hallmark,
  gs_C2_kegg = gs_C2_kegg,
  gs_C2_reactome = gs_C2_reactome,
  gs_C5_GOBP = gs_C5_GOBP,
  gs_C5_GOCC = gs_C5_GOCC, 
  gs_C5_GOMF = gs_C5_GOMF
)

pdf(file = file.path(deg_dir, "MES_vs_ADR_fGSEA_and_dotPLOT2.pdf"), width = 10, height = 10)
for(gene_set_name in names(gene_set_list)) {
  print(gene_set_name)
  gene_set <- gene_set_list[[gene_set_name]]
  
  fgseaRes <- fgsea(pathways = gene_set[["genesets"]], 
                    stats    = de_genes_ranked,
                    minSize  = 15,
                    maxSize  = 500)
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  p <- plotGseaTable(gene_set[["genesets"]][topPathways],
                     de_genes_ranked,
                     fgseaRes, 
                     gseaParam = 0.5) + 
    labs(title = gene_set_name)
  p[["theme"]][["plot.title"]] <- NULL
  print(p)
  
  # Dot plot
  fgseaRes_dot <- fgseaRes %>% 
    select(pathway, NES, padj) %>%
    filter(pathway  %in% topPathways) %>%
    arrange(NES) %>%
    mutate(pathway=factor(pathway, levels = pathway)) %>%
    mutate(padj_log2 = -log2(padj))
  
  p <- ggplot(fgseaRes_dot) +
    aes(x = NES, y = pathway, colour = NES, size = padj_log2) +
    geom_point(shape = "circle") +
    scale_color_distiller(palette = "RdBu", direction = -1) +
    scale_size(range = c(4, 10)) +
    theme_minimal() +
    labs(title = gene_set_name,
         size = "-log2p_adj")
  print(p)
}
dev.off()


# Create a volcano plot
# Genes to highlight
genes_to_highlight <- c("VIM",
                        "YAP1",
                        "WWTR1",
                        "JUN",
                        "FOSL1",
                        "FOSL2",
                        "PHOX2B",
                        "HAND2",
                        "GATA3")
volcano_plot <- plotVolcano(dds_results_obj = deg_results$results_all, 
                                                    genes_of_interest = genes_to_highlight, 
                                                    plot_title = "MES vs ADRN genes")
ggsave(filename = paste0(deg_dir, "MES_vs_ADRN_volcano_plot.png"), 
       plot = volcano_plot,
       width = 20, height = 20, units = "cm")






# # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # EOF         # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # 
# perform ssGSEA on Thirant data
thir_data <- read.table(file = "~/workspace/neuroblastoma/resources/Thiriant_data.csv", 
                        sep = ";",
                        header = TRUE)

thir_data <- split(thir_data %>% select(group, gene_name), thir_data$group)
thir_data <- lapply(thir_data, function(x) return(x$gene_name))

fgseaRes <- fgsea(pathways = thir_data["OX-PHOS"], 
                  stats    = de_genes_ranked,
                  minSize  = 15,
                  maxSize  = 500)
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
p <- plotGseaTable(gs_C5_GOBP[["genesets"]][topPathways],
                   de_genes_ranked,
                   fgseaRes, 
                   gseaParam=0.5) + labs(title = "C5_GOBP")
p[["theme"]][["plot.title"]] <- NULL
p


###### DEV ZONE
# Plot terms based on terms defined from Thirant supplementary data 1:
thir_data <- read.table(file = "~/workspace/neuroblastoma/resources/Thirant_supplementary_data.tsv",
                       sep = "\t",
                       header = TRUE)
thir_data$group_cluster <- paste0(thir_data$group, "_", thir_data$cluster)

thir_data_group <- base::split(thir_data, thir_data$group)
thir_data_group <- lapply(thir_data_group, function(x){pull(x, gene_name)})

thir_data_cluster <- base::split(thir_data, thir_data$group_cluster)
thir_data_cluster <- lapply(thir_data_cluster, function(x){pull(x, gene_name)})

GSEA_thir_data <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                        genesets = thir_data_group, 
                        test="hypergeometric", 
                        background=nrow(rnaseq_dds_filt))
GSEA_thir_data_plot <- hypeR::hyp_dots(GSEA_thir_data, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "Thirant sup. data: MES vs ADR") +
  theme_bw()

GSEA_thir_data <- hypeR::hypeR(signature = deg_results$results_signif$gene_symbol, 
                               genesets = thir_data_cluster, 
                               test = "hypergeometric", 
                               background = nrow(rnaseq_dds_filt))
GSEA_thir_data_plot <- hypeR::hyp_dots(GSEA_thir_data, 
                                       merge = TRUE, 
                                       fdr = 0.05, 
                                       top = 20, 
                                       abrv = 70, 
                                       val = "fdr", 
                                       title = "Thirant sup. data: MES vs ADR") +
  theme_bw()

#####
thir_data <- read.table(file = "~/workspace/neuroblastoma/resources/Thirant_supplementary_data.tsv",
                        sep = "\t",
                        header = TRUE)
thir_data <- read.table(file = "~/workspace/neuroblastoma/resources/Thiriant_data.csv",
                        sep = ";",
                        header = TRUE)

thir_data$group_cluster <- paste0(thir_data$group, "_", thir_data$cluster)

thir_data_group <- base::split(thir_data, thir_data$group)
thir_data_group <- lapply(thir_data_group, function(x){pull(x, gene_name)})

thir_data_cluster <- base::split(thir_data, thir_data$group_cluster)
thir_data_cluster <- lapply(thir_data_cluster, function(x){pull(x, gene_name)})


vsd_counts_matrix <- assay(vsd)
row.names(vsd_counts_matrix) <- annotationData[, 'gene_symbol'][match(row.names(vsd_counts_matrix), annotationData[, 'ensembl_id'])]

ssgsea_mes_adr_thir_data<- GSVA::gsva(vsd_counts_matrix,
                                      thir_data_group,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      #ssgsea.norm=TRUE, 
                                      verbose=TRUE, 
                                      parallel.sz=10,
                                      )
ssgsea_mes_adr_thir_data <- as.data.frame(ssgsea_mes_adr_thir_data)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_thir_data,
                                                        scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = FALSE,
                                                        cluster_cols = TRUE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)


ssgsea_mes_adr_thir_data_cluster<- GSVA::gsva(vsd_counts_matrix,
                                      thir_data_cluster,
                                      method=c("ssgsea"),
                                      min.sz=1, max.sz=Inf, 
                                      #ssgsea.norm=TRUE, 
                                      verbose=TRUE, 
                                      parallel.sz=10,
)
ssgsea_mes_adr_thir_data_cluster <- as.data.frame(ssgsea_mes_adr_thir_data_cluster)
ssgsea_mes_adr_ncc_noradr_heatmap <- pheatmap::pheatmap(ssgsea_mes_adr_thir_data_cluster,
                                                        scale = "row",
                                                        #annotation_col = heatmap_col_annot_DTC,
                                                        cluster_rows = FALSE,
                                                        cluster_cols = TRUE,
                                                        color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                                                        show_colnames = TRUE)

#############




# Produce heat map - all DE genes in ADR and MES lines 
metadata_heatmap <- as.data.frame(colData(rnaseq_dds_filt))
heatmap_counts <- SummarizedExperiment::assay(vsd)
heatmap_counts_deg <- heatmap_counts[rownames(heatmap_counts) %in% deg_results$results_signif$ensembl_id, ]
annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_type, cell_line) %>% 
  dplyr::arrange(cell_type, cell_line)
heatmap_counts_deg_ord <- heatmap_counts_deg[, match(rownames(annotation_col), colnames(heatmap_counts_deg))]
color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
#color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
ann_colors = list(
  cell_line = c(CLB_Ma = "#005f73", STA_NB_10 = "#0a9396", STA_NB_8 = "#DD3344", SK_N_SH = "#FF9F1C"),
  cell_line_id = c(ADR = "#E9D8A6", MES = "#D9BE6D")
)
heatmap <- pheatmap::pheatmap(heatmap_counts_deg_ord,
                                            main = "Heatmap of signif. DEG",
                                            scale = "row",
                                            annotation_col = annotation_col,
                                            annotation_colors = ann_colors,
                                            show_colnames = FALSE,
                                            show_rownames = FALSE,
                                            cluster_cols = FALSE,
                                            color = color.scheme,
                                            fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
ggsave(filename = paste0(deg_dir, "ADRN_vs_MES_heatmap.png"), 
       plot=heatmap,
       width = 20, height = 20, units = "cm")





### Extract counts for diffTF ###
Counts_for_diffTF <- rnaseq_dds_filt@assays@data@listData[["counts"]]
row.names(Counts_for_diffTF) <- row.names(rnaseq_dds_filt)
colnames(Counts_for_diffTF) <- c( gsub("_RNA_S.*", replacement = "", colnames(Counts_for_diffTF)))
#colnames(Counts_for_diffTF) <- c("ENSEMBL", gsub("_RNA_S.*", replacement = "", colnames(Counts_for_diffTF)))
write.table(Counts_for_diffTF, file = file.path(deg_dir, "RNAseq.tsv"), 
            sep = "\t", 
            row.names = TRUE,
            col.names = TRUE,
            quote = FALSE)









### xcore ####
#promoters_f5_core <- xcoredata::promoters_f5_core()

#taken from the previous steps
dds_counts <- counts(rnaseq_dds_filt, normalized=FALSE)

#subset the pattern to create Design table
cond <- str_extract(string = colnames(dds_counts),
                    pattern = "_(A|M)_", 
                    group = 1)

design <- data.frame(row.names = colnames(dds_counts))
design$A <- ifelse(test = cond == "A", 
                   yes = 1,
                   no = 0)
design$M <- ifelse(test = cond == "M", 
                   yes = 1,
                   no = 0)
design <- as.matrix(design)

# load F5 data and symbol2fantom
promoters_f5_core <- xcoredata::promoters_f5_core()
remap_promoters_f5 <- xcoredata::remap_promoters_f5()

eh <- ExperimentHub::ExperimentHub()
symbol2fantom <- eh[["EH7700"]]

# using annotationData from previous step
annotationData

# replace ENSEMBL IDs with symbols and remove duplicates
row.names(dds_counts) <-  annotationData$gene_symbol[match(row.names(dds_counts) , annotationData$ensembl_id)]
dds_counts <- dds_counts[-which(duplicated(row.names(dds_counts))),]
counts_rna_seq_fantom <- translateCounts(dds_counts, dict = symbol2fantom)

# main xcore part
mae_rna_seq <- prepareCountsForRegression(
  counts = counts_rna_seq_fantom,
  design = design,
  base_lvl = "A"
)
mae_rna_seq <- addSignatures(mae_rna_seq, 
                             remap = remap_promoters_f5)
mae_rna_seq <- filterSignatures(mae_rna_seq, min = 0.05, max = 0.95)

# register parralel backend
doMC::registerDoMC(cores = 6L)
# set seed
set.seed(314159265)

res_rna_seq <- modelGeneExpression(
  mae = mae_rna_seq,
  xnames = "remap",
  nfolds = 6
  )

#find what factors have pval < 0.05
list_of_factors <- list()
list_of_factors_bool <- data.frame(row.names = res_rna_seq[["results"]][["remap"]][["name"]])

for(sample_name in names(res_rna_seq[["regression_models"]][["remap"]])){
  list_of_factors[[sample_name]] <- data.frame(res_rna_seq[["pvalues"]][["remap"]][[sample_name]]) %>% 
                                    dplyr::select("pval") %>% 
                                    dplyr::filter(pval < 0.05) %>% 
                                    row.names()
  list_of_factors_bool[,sample_name] <- res_rna_seq[["pvalues"]][["remap"]][[sample_name]][["pval"]] < 0.05
}

#create a bool vector
list_of_factors_bool$n_pval_filt <- (apply(list_of_factors_bool, 1, sum))
at_least_one_motif <- list_of_factors_bool %>% dplyr::filter(n_pval_filt > 0) %>% row.names()
unique_motifs <- Reduce(intersect, list_of_factors)

top_signatures  <- res_rna_seq$results$remap %>% dplyr::filter(name %in% unique_motifs)
pheatmap::pheatmap(
  mat = top_signatures[, "M"],
  labels_row = top_signatures$name,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(35),
  breaks = seq(from = -0.2, to = 0.2, length.out = 36),
  main = "ReMap2020 molecular signatures activity"
)

top_signatures2  <- res_rna_seq$results$remap %>% dplyr::filter(name %in% at_least_one_motif)
pheatmap::pheatmap(
  mat = top_signatures2[, "M"][1:100],
  labels_row = top_signatures2$name[1:100],
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(35),
  breaks = seq(from = -0.2, to = 0.2, length.out = 36),
  main = "ReMap2020 molecular signatures activity"
)
