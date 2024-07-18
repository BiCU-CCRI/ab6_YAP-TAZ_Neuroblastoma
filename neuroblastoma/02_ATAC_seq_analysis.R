###
# Title: ATAC-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
# 

# Set up the environment ####
deg_dir <- "~/workspace/neuroblastoma/results/ATAC-seq/"

## Loading libraries ####
library(dplyr)
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
             "extract_results_DDS",
             "meanExprsPerGroup",
             "extract_results_DDS_HIC",
             "plotVolcano",
             "get_the_poissoon_p_val")
# Bug - JASPAR2022 cannot be pulled from the web automatically, so it need to be loaded manually
download.file(url = "https://jaspar2022.genereg.net/download/database/JASPAR2022.sqlite",
              destfile = "/home/rstudio/.cache/R/JASPAR2022.sqlite")
JASPAR2022 <-  "/home/rstudio/.cache/R/JASPAR2022.sqlite"

# Analysis part ####
## Assembling the metadata to a DESeq2 object ####

# Assembling the metadata to a DESeq2 object
# we use only filtered consensus peaks. And combined annotation from genecode, homer and reg.
tmp_dir <- "~/workspace/neuroblastoma/data/ATACseq/"
path_to_the_ATAC_counts_data <- file.path(tmp_dir, "all_filtered.RDS")
path_annnotation_ATAC_data <- file.path(tmp_dir, "consensus_regions_annotation.RDS")
ATAC_counts_data <- readRDS(path_to_the_ATAC_counts_data)
ATAC_annotation_data <- readRDS(path_annnotation_ATAC_data)

#subset the ATAC_annotation_data - keep the rows that are present in the ATAC_counts_data
ATAC_annotation_data <- dplyr::filter(ATAC_annotation_data, peak_id %in% rownames(ATAC_counts_data))
# prepare metadata table form the names of the samples
ATAC_metadata_df <- data.frame(matrix(nrow = dim(ATAC_counts_data)[2]))
ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
                                  sample_names = colnames(ATAC_counts_data)
)
ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
                                  cell_line = str_extract(pattern = "^(.*?)_", 
                                                          string = colnames(ATAC_counts_data)
                                  )
)
ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
                                  cell_line = gsub(pattern = "_", 
                                                   replacement = "", 
                                                   x = ATAC_metadata_df$cell_line
                                  )
)
ATAC_metadata_df <- ATAC_metadata_df[, -1]
ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
                                  phenotype = str_extract(pattern = "_._", 
                                                          string = ATAC_metadata_df$sample_names
                                  )
)
ATAC_metadata_df <- dplyr::mutate(ATAC_metadata_df, 
                                  phenotype = gsub(pattern = "_", 
                                                   replacement = "", 
                                                   x = ATAC_metadata_df$phenotype)
)
# Assembling the DESeq2 object. 
desing_formula <- as.formula("~ cell_line + phenotype")

# The Warning "some variables in design formula are characters, converting to factors" is normal.
ATAC_dds <- DESeq2::DESeqDataSetFromMatrix(countData = ATAC_counts_data, 
                                           colData = ATAC_metadata_df, 
                                           design = desing_formula
)
#Adding ranges info
rowRanges(ATAC_dds) <- GenomicRanges::GRanges(seqnames = ATAC_annotation_data$gencode_chr, 
                                              ranges = IRanges::IRanges(start = ATAC_annotation_data$gencode_start, 
                                                                        end = ATAC_annotation_data$gencode_end
                                              )
)
# Adding metadata
mcols(ATAC_dds) <- cbind(mcols(ATAC_dds), ATAC_annotation_data)        

# plot full ATAC-seq PCA
ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)
vsd <- vst(ATAC_dds, blind = F)

# Generate PCA
pca_dar <- generatePCA(transf_object = vsd, 
                       cond_interest_varPart = c("phenotype", "cell_line"), 
                       color_variable = "phenotype", 
                       shape_variable = "cell_line",
                       ntop_genes = 1000) +
  ggtitle("Original dataset") +
  scale_color_manual(values = c("#DD3344",  "#FF9F1C"))
ggsave(
  filename = file.path(deg_dir, "PCA_MES_vs_ADRN_DARs_all_samples.png"),
  plot = pca_dar,
  width = 20, height = 20, units = "cm"
)

# We remove NB6 samples, as they are mixed cultures of ADRN and MES
ATAC_dds <- ATAC_dds[, !(colnames(ATAC_dds) %in% c("NB6_A_ATAC", "NB6_M_ATAC"))]
ATAC_dds$cell_line <- droplevels(ATAC_dds$cell_line)

# # Load the H3K27Ac data and converting the coordinates to hg38
# # Important REMARK:
# # Importing BED files with H3K27Ac (labelling super enhancers)  that are form Groningen paper. 
# # Weirdly enough - They claim that they managed to identify 286 MES specific and 276 ADRN-specific sEnh. 
# # But when I try to filter the table that they provide using FDR <0.1 I get 254 and 255 sEnh correspondingly. 
# 
# # BED_colnames <- c("Chr", "Start", "End")
# # Boevaetal_ADRN_H3K27Ac_ChiP <- read.table(file = "./neuroblastoma/data_public/data/Boevaetal_ADRN_H3K27Ac_ChiP.bed", header = F, sep = "\t" )[,1:3]
# # Boevaetal_MES_H3K27Ac_ChiP <- read.table(file = "./neuroblastoma/data_public/data/Boevaetal_MES_H3K27Ac_ChiP.bed", header = F, sep = "\t")[,1:3]
# 
# Groningenetal_H3K27Ac <- read_xlsx(xlsxFile = "./neuroblastoma/data_public/data/Groningenetal_superEnh.xlsx",startRow = 3)
# Groningenetal_H3K27Ac <- dplyr::filter(Groningenetal_H3K27Ac,  limma_fdr < 0.1 & limma_pval < 0.05)
# Groningenetal_ADRN_H3K27Ac <- Groningenetal_H3K27Ac[str_detect(string = Groningenetal_H3K27Ac$sample_ids, pattern = "ADRN") & 
#                                                       !(str_detect(string = Groningenetal_H3K27Ac$sample_ids, pattern = "MES")), 1:4]
# Groningenetal_MES_H3K27Ac <- Groningenetal_H3K27Ac[str_detect(string = Groningenetal_H3K27Ac$sample_ids, pattern = "MES") & 
#                                                       !(str_detect(string = Groningenetal_H3K27Ac$sample_ids, pattern = "ADRN")), 1:4]
# 
# # Convert hg19 coordinates to hg38
# # Specify coordinates to liftover
# Groningenetal_ADRN_H3K27Ac <- GRanges(Groningenetal_ADRN_H3K27Ac)
# Groningenetal_MES_H3K27Ac <- GRanges(Groningenetal_MES_H3K27Ac)
# # Import the chain file
# chainObject <- import.chain("./neuroblastoma/resources/hg19ToHg38.over.chain")
# # Run the liftOver
# Groningenetal_ADRN_H3K27Ac <- GRanges(as.data.frame(liftOver(Groningenetal_ADRN_H3K27Ac, chainObject)))
# Groningenetal_MES_H3K27Ac <- GRanges(as.data.frame(liftOver(Groningenetal_MES_H3K27Ac, chainObject)))
# 
# # Now, as we have super enhancer regions, we can filter the regions in our ATAC-seq data and add this metadata to our atac-seq dataset 
# ADRN_peak_ids <- subsetByOverlaps(rowRanges(ATAC_dds), 
#                                   Groningenetal_ADRN_H3K27Ac
#                                   )$peak_id
# MES_peak_ids <- subsetByOverlaps(rowRanges(ATAC_dds), 
#                                  Groningenetal_MES_H3K27Ac
#                                  )$peak_id
# mcols(ATAC_dds)$sEnh_type <- "ND"
# mcols(ATAC_dds)$sEnh_type[mcols(ATAC_dds)$peak_id %in% ADRN_peak_ids] <- "ADRN_sh"
# mcols(ATAC_dds)$sEnh_type[mcols(ATAC_dds)$peak_id %in% MES_peak_ids] <- "MES_sh"
# 
# # Mark peaks that belong to genes that are known to be differentially expressed in ADRN and MES
# ADRN_MES_genes <- read_xlsx(xlsxFile = "./neuroblastoma/resources/41588_2017_BFng3899_MOESM3_ESM.xlsx",
#                             startRow = 2
#                             )
# colnames(ADRN_MES_genes) <- c("Gene", "category")
# mcols(ATAC_dds)$gene_category_GROEN <- "ND"
# mcols(ATAC_dds)$gene_category_GROEN <- case_when(
#   mcols(ATAC_dds)$homer_Gene.Name %in% ADRN_MES_genes[ADRN_MES_genes$category == "ADRN", 1] ~ "g_ADRN",
#   mcols(ATAC_dds)$homer_Gene.Name %in% ADRN_MES_genes[ADRN_MES_genes$category == "MES", 1] ~ "g_MES",
#   TRUE ~ mcols(ATAC_dds)$gene_category_GROEN
# )
# 
# Loading OUR own rna-seq data and marking ADRN and MES-specific genes

RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
RNA_SEQ_data <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES")) %>%
  select(gene_symbol, Term)
mcols(ATAC_dds)$gene_category_our_RNAseq <- "ND"
mcols(ATAC_dds)$gene_category_our_RNAseq <- case_when(
  mcols(ATAC_dds)$homer_Gene.Name %in% RNA_SEQ_data[RNA_SEQ_data$Term == "ADRN", 1] ~ "g_ADRN",
  mcols(ATAC_dds)$homer_Gene.Name %in% RNA_SEQ_data[RNA_SEQ_data$Term == "MES", 1] ~ "g_MES",
  TRUE ~ mcols(ATAC_dds)$gene_category_our_RNAseq
)

# Add metadata - marks peaks that are around 3kb from TSS
mcols(ATAC_dds)$is_promoter_3kb <- FALSE
mcols(ATAC_dds)$is_promoter_3kb[abs(mcols(ATAC_dds)$homer_Distance.to.TSS) <= 3000] <- TRUE

# # Add metadata - if peaks are in the area of YAP TAZ
# Groningenetal_WWTR <- read.table(file = "~/workspace/neuroblastoma/data_public/data/Groningenetal_WWTR-YAP.bed", sep = "\t")
# colnames(Groningenetal_WWTR) <- c("chr", "start", "end")
# Groningenetal_WWTR <- unlist(liftOver(GRanges(Groningenetal_WWTR), chainObject))
# Boevaetal_WWTR <- read.table(file = "~/workspace/neuroblastoma/data_public/data/Boevaetal_WWTR-YAP.bed", sep = "\t")
# colnames(Boevaetal_WWTR) <- c("chr", "start", "end")
# Boevaetal_WWTR <- unlist(liftOver(GRanges(Boevaetal_WWTR), chainObject))
# 
# # Combine these ranges into one object and subset it 
# TAZ_ranges <- reduce(unlist(GRangesList(Groningenetal_WWTR , Boevaetal_WWTR)), 
#                      drop.empty.ranges = T, 
#                      min.gapwidth = 1
#                      )[1:4,]
# YAP_ranges <- reduce(unlist(GRangesList(Groningenetal_WWTR,Boevaetal_WWTR)), 
#                      drop.empty.ranges = T, 
#                      min.gapwidth = 1
#                      )[5:7,]
# # Overlap the peaks location with these ranges
# YAP_ids <- subsetByOverlaps(rowRanges(ATAC_dds), YAP_ranges)$peak_id
# TAZ_ids <- subsetByOverlaps(rowRanges(ATAC_dds), TAZ_ranges)$peak_id
# mcols(ATAC_dds)$YAP_TAZ <- "ND"
# mcols(ATAC_dds)$YAP_TAZ[mcols(ATAC_dds)$peak_id %in% YAP_ids] <- "YAP"
# mcols(ATAC_dds)$YAP_TAZ[mcols(ATAC_dds)$peak_id %in% TAZ_ids] <- "TAZ"
# 
# # Save the object
# saveRDS(ATAC_dds, file = "~/workspace/neuroblastoma/RDSs/ATAC_dds.RDS")


#### General analysis fo the ATAC-seq ####
# Metadata is added, now we can process the ATAC-seq data set and do the DE +
# variance stabilization
#ATAC_dds <- ATAC_dds_full
ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)
vsd <- vst(ATAC_dds, blind = F)

# Generate PCA
pca_dar <- generatePCA(transf_object = vsd, 
                       cond_interest_varPart = c("phenotype", "cell_line"), 
                       color_variable = "phenotype", 
                       shape_variable = "cell_line",
                       ntop_genes = 1000) +
  ggtitle("Original dataset") +
  scale_color_manual(values = c("#DD3344",  "#FF9F1C"))
ggsave(
  filename = file.path(deg_dir, "PCA_MES_vs_ADRN_DARs.png"),
  plot = pca_dar,
  width = 20, height = 20, units = "cm"
)
# Looks good - we can see that the main PCA/separation is happening because of the phenotype (ADRN or MES)

# Extracting results
resultsNames(ATAC_dds)
coeff_name <- "phenotype_M_vs_A"
cond_numerator <-  "A"
cond_denominator <-  "M"
cond_variable <- "phenotype"
padj_cutoff = 0.05
log2FC_cutoff = 1

ATAC_dds_results <- extract_results_DDS(dds_object = ATAC_dds,
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

openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_all, sheet = "results_all")

openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir, "DAR_results_M_VS_A.xlsx"), overwrite = T)

# making volcanoplot
genes_up <- c("VIM", "FOSL2", "FOSL1", "YAP1", "JUN", "WWTR1")
genes_down <- c("PHOX2B", "HAND2", "GATA3")
p1 <- plotVolcano(dds_results_obj = ATAC_dds_results$results_all,
                  genes_of_interest = c("VIM", "FOSL2", "FOSL1", "YAP1", "JUN", "WWTR1", "PHOX2B", "HAND2", "GATA3"))
p1

# Making Heatmap
metadata_heatmap <- as.data.frame(colData(ATAC_dds))
heatmap_counts <- SummarizedExperiment::assay(vsd)
rownames(heatmap_counts) <- vsd@rowRanges$peak_id

# Subset different peaks that belong to different categories
ATAC_signif <- ATAC_dds_results$results_signif
heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id, ]
#heatmap_counts_ADRN <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$sEnh_type == "ADRN_sh"], ]
#heatmap_counts_MES <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$sEnh_type == "MES_sh"], ]
#heatmap_counts_YAP <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$YAP_TAZ == "YAP"],]
#heatmap_counts_TAZ <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$YAP_TAZ == "TAZ"],]
heatmap_counts_MES_spec_genes <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$gene_category_our_RNAseq == "g_MES"],]
heatmap_counts_ADRN_spec_genes <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$gene_category_our_RNAseq == "g_ADRN"],]
  
  
annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, phenotype) %>% 
  dplyr::arrange(phenotype, cell_line)
heatmap_counts<- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]
# heatmap_counts_ADRN <- heatmap_counts_ADRN[, match(rownames(annotation_col), colnames(heatmap_counts_ADRN))]
# heatmap_counts_MES <- heatmap_counts_MES[, match(rownames(annotation_col), colnames(heatmap_counts_MES))]
# heatmap_counts_YAP <- heatmap_counts_YAP[, match(rownames(annotation_col), colnames(heatmap_counts_YAP))]
# heatmap_counts_TAZ <- heatmap_counts_TAZ[, match(rownames(annotation_col), colnames(heatmap_counts_TAZ))]
heatmap_counts_MES_spec_genes <- heatmap_counts_MES_spec_genes[, match(rownames(annotation_col), colnames(heatmap_counts_MES_spec_genes))]
heatmap_counts_ADRN_spec_genes <- heatmap_counts_ADRN_spec_genes[, match(rownames(annotation_col), colnames(heatmap_counts_ADRN_spec_genes))]
  
ensembl2symbol_annot <- ATAC_dds_results$results_all %>%
  dplyr::select(peak_id, gencode_gene_name)

#Alternative color schemes
#color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012", "#9B2226")
#color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
#color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
color.scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(50)
ann_colors = list(
  cell_line = c(CLBM = "#005f73", NB10 = "#0a9396", NB8 = "#DD3344", SH = "#FF9F1C"),
  cell_line_id = c(A = "#E9D8A6", M = "#D9BE6D")
)

##### exporting clustered peaks to a separate fiels
# This is for heatmap-denstify plots 
hclust_matrix <- heatmap_counts %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()
gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
plot(gene_hclust, labels = FALSE)
clustered_peaks <- cutree(gene_hclust, k = 2)
df_clustered_peaks <- data_frame(names(clustered_peaks), clustered_peaks)
colnames(df_clustered_peaks) <- c("peak_id", "cluster")
df_clustered_peaks <- merge(x = ATAC_dds_results$results_signif, 
      y = df_clustered_peaks,
      by = "peak_id")
df_clustered_peaks <- df_clustered_peaks %>% dplyr::select(gencode_chr, gencode_start, gencode_end, peak_id, cluster)
write.table(
  df_clustered_peaks %>% 
    dplyr::filter(cluster == 1) %>%
    dplyr::select(-cluster),
  file = "~/workspace/neuroblastoma/tmp/tmp_bed_cluster_1.bed",
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
write.table(
  df_clustered_peaks %>% 
    dplyr::filter(cluster == 2) %>%
    dplyr::select(-cluster),
  file = "~/workspace/neuroblastoma/tmp/tmp_bed_cluster_2.bed",
  sep = "\t",
  quote = FALSE, 
  row.names = FALSE, 
  col.names = FALSE
)
###########

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "Heatmap of signif. DAR ADRN vs MES.",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

#ATAC_dds_results$results_signif[ATAC_dds_results$results_signif$log2FoldChange > 0,] 

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_full_ATAC.png"), 
       plot = heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_full_ATAC.eps"), 
       plot = heatmap,
       width = 18, height = 20, units = "cm")

# # Now we can check how the peaks accessibility in MES/ADRN-specific regions identified in previous studies
# heatmap <- pheatmap::pheatmap(heatmap_counts_ADRN,
#                               main = "Heatmap of signif. DAR ADRN vs MES. ADRN-specific",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) 
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN_enh_ATAC.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN_enh_ATAC.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# 
# heatmap <- pheatmap::pheatmap(heatmap_counts_MES,
#                               main = "Heatmap of signif. DAR ADRN vs MES. MES-specific",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               #annotation_row = row_annot_symbols,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_enh_ATAC.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_enh_ATAC.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# 
# heatmap <- pheatmap::pheatmap(heatmap_counts_YAP,
#                               main = "Heatmap of signif. DAR ADRN vs MES. YAP",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               #annotation_row = row_annot_symbols,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_YAP.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_YAP.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# 
# heatmap <- pheatmap::pheatmap(heatmap_counts_TAZ,
#                               main = "Heatmap of signif. DAR ADRN vs MES. TAZ",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               #annotation_row = row_annot_symbols,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_TAZ.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_TAZ.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")


heatmap <- pheatmap::pheatmap(heatmap_counts_MES_spec_genes,
                              main = "Heatmap of signif. DAR ADRN vs MES. MES-specific genes only (defined by RNA-seq)",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES-RNAseq-specific.png"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES-RNAseq-specific.eps"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")


heatmap <- pheatmap::pheatmap(heatmap_counts_ADRN_spec_genes,
                              main = "Heatmap of signif. DAR ADRN vs MES. ADRN-specific genes only (defined by RNA-seq)",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN-RNAseq-specific.png"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN-RNAseq-specific.eps"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")


## Genomic features enrichment lolipop plot ######
# use HOMER's output
# Separate positive - M negative - A
ATAC_dds_results$results_signif %>%
  filter(log2FoldChange > 0) %>%
  select(gencode_chr, gencode_start, gencode_end) %>%
  write.table(., file = paste0(deg_dir, "/bed_diff_MES.bed"),
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE,
              col.names = FALSE)

ATAC_dds_results$results_signif %>%
  filter(log2FoldChange < 0) %>%
  select(gencode_chr, gencode_start, gencode_end) %>%
  write.table(., file = paste0(deg_dir, "/bed_diff_ADR.bed"),
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE,
              col.names = FALSE)

# run homer annotation
system("bash ~/workspace/neuroblastoma/02_HOMER_diff_peaks_annotation.sh")

ATAC_genome_ADR <- read.delim2("~/workspace/neuroblastoma/results/ATAC-seq/ATAC_genome_states_dataframe_ADR.csv", 
                               sep="\t", 
                               header=T)[-c(1:14),]
ATAC_genome_MES <- read.delim2("~/workspace/neuroblastoma/results/ATAC-seq/ATAC_genome_states_dataframe_MES.csv",
                               sep="\t",
                               header=T)[-c(1:14),]

# calculate the poisson probability


ATAC_genome_ADR$p.val <- unlist(
  purrr::pmap(ATAC_genome_ADR %>% 
                dplyr::select(Number.of.peaks, Log2.Ratio..obs.exp.),
              get_the_poissoon_p_val) 
)
ATAC_genome_ADR$p.adj <- p.adjust(ATAC_genome_ADR$p.val, method = "BH")
ATAC_genome_ADR$p.adj <- round(ATAC_genome_ADR$p.adj, 6)
ATAC_genome_ADR$p.adj[ATAC_genome_ADR$p.adj == 0] <- 0.000001
ATAC_genome_ADR$negLog2p.adj <- -log2(ATAC_genome_ADR$p.adj)
ATAC_genome_ADR <- ATAC_genome_ADR[,c(1,4,8)]

ATAC_genome_MES$p.val <- unlist(
  purrr::pmap(ATAC_genome_MES %>% 
                dplyr::select(Number.of.peaks, Log2.Ratio..obs.exp.),
              get_the_poissoon_p_val) 
)
ATAC_genome_MES$p.adj <- p.adjust(ATAC_genome_MES$p.val, method = "BH")
ATAC_genome_MES$p.adj <- round(ATAC_genome_MES$p.adj, 6) 
ATAC_genome_MES$p.adj[ATAC_genome_MES$p.adj == 0] <- 0.000001
ATAC_genome_MES$negLog2p.adj <- -log2(ATAC_genome_MES$p.adj)
ATAC_genome_MES <- ATAC_genome_MES[,c(1,4,8)]

ATAC_genome <- merge.data.frame(ATAC_genome_ADR, ATAC_genome_MES, by = "Annotation", suffixes = c(".ADR", ".MES"))
ATAC_genome$Log2.Ratio..obs.exp..ADR <- as.numeric(ATAC_genome$Log2.Ratio..obs.exp..ADR)
ATAC_genome$Log2.Ratio..obs.exp..MES <- as.numeric(ATAC_genome$Log2.Ratio..obs.exp..MES)

df_long <- ATAC_genome %>%
  tidyr::pivot_longer(
    cols = starts_with("Log2.Ratio") | starts_with("negLog2p.adj"),
    names_to = c(".value", "Type"),
    names_pattern = "(Log2.Ratio..obs.exp..|negLog2p.adj.)(ADR|MES)"
  )

pdf(file = file.path(deg_dir, "MES_and_ADR_lolipop_enrichemtn_plot.pdf"), width = 8, height = 8)
df_long %>%
  filter(Annotation %in% c("3UTR", "5UTR", "Exon", "Intergenic", "Intron", "Promoter")) %>%
  ggplot() +
  geom_bar(position = position_dodge(0.5), 
           width = 0.1, 
           aes(y = Annotation, 
               fill = Type, 
               weight = Log2.Ratio..obs.exp..)) +
  scale_fill_manual(values =  c("navy", "firebrick3")) +
  #  scale_fill_hue(direction = 1) +
  geom_point(aes(y = Annotation, 
                 x = Log2.Ratio..obs.exp.., 
                 color = Type,
                 size = negLog2p.adj.), 
             shape = 21, 
             position = position_dodge2(0.5),
             fill = "white") +
  scale_color_manual(values =  c("navy", "firebrick3")) +
  #  scale_size_continuous(range = c(2, 10)) +
  #  scale_color_gradient(low = "blue", high = "red", limits = c(0,20)) +
  theme_minimal() +
  labs(color = "-Log2p.adj.", 
       fill = "Type",
       size =  "-Log2p.adj.",
       y = "Annotation",
       x = "Log2 Ratio (obs/exp)",
       title = "Enrichment of genomic features in ATAC-seq peaks for ADRN and MES")

dev.off()






#make a small back up 
ATAC_dds_full <- ATAC_dds

#### now the same analysis, but for the promoter regions - ± 3000bp around TSSs ####
ATAC_dds <- ATAC_dds[mcols(ATAC_dds)$is_promoter_3kb,]
ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)

vsd <- vst(ATAC_dds, blind = F)

pca_dar <- generatePCA(transf_object = vsd, 
                       cond_interest_varPart = c("phenotype", "cell_line"), 
                       color_variable = "phenotype", 
                       shape_variable = "cell_line",
                       ntop_genes = 1000) +
  ggtitle("3kb Promoters") +
  scale_color_manual(values = c("#DD3344",  "#FF9F1C"))
plot(pca_dar)
#Ok, looks good - we can see that the main PCA/separation is happening because of the phenotype (ADRN or MES)

# extracting results

resultsNames(ATAC_dds)
coeff_name <- 'phenotype_M_vs_A'
cond_numerator <-  "A"
cond_denominator <-  "M"
cond_variable <-  "phenotype"
padj_cutoff = 0.05
log2FC_cutoff = 1

ATAC_dds_results <- extract_results_DDS(dds_object = ATAC_dds,
                                        coeff_name = coeff_name,
                                        cond_numerator = cond_numerator,
                                        cond_denominator = cond_denominator,
                                        cond_variable = cond_variable,
                                        padj_cutoff = padj_cutoff,
                                        log2FC_cutoff = log2FC_cutoff)

#write the result to the xlsx file
XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")

openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_all, sheet = "results_all")

openxlsx::saveWorkbook(XLSX_OUT, file.path(deg_dir, "DAR_results_3kb_promoters_M_VS_A.xlsx"), overwrite = T)

#Making Heatmap
metadata_heatmap <- as.data.frame(colData(ATAC_dds))
heatmap_counts <- SummarizedExperiment::assay(vsd)
rownames(heatmap_counts) <- vsd@rowRanges$PeakId

heatmap_counts <- SummarizedExperiment::assay(vsd) 
rownames(heatmap_counts) <- vsd@rowRanges$peak_id
# removed experiment batch effects! ? use experimen+cell_line removed???

heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_dds_results$results_signif$peak_id,]
# heatmap_counts_ADRN <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_dds_results$results_signif$peak_id[ATAC_dds_results$results_signif$sEnh_type == "ADRN_sh"], ]
# heatmap_counts_MES <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_dds_results$results_signif$peak_id[ATAC_dds_results$results_signif$sEnh_type == "MES_sh"], ]
heatmap_counts_MES_spec_genes <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$gene_category_our_RNAseq == "g_MES"],]
heatmap_counts_ADRN_spec_genes <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_signif$peak_id[ATAC_signif$gene_category_our_RNAseq == "g_ADRN"],]

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, phenotype) %>% 
  dplyr::arrange( phenotype, cell_line)
heatmap_counts<- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]
# heatmap_counts_ADRN<- heatmap_counts_ADRN[, match(rownames(annotation_col), colnames(heatmap_counts_ADRN))]
# heatmap_counts_MES<- heatmap_counts_MES[, match(rownames(annotation_col), colnames(heatmap_counts_MES))]
heatmap_counts_MES_spec_genes <- heatmap_counts_MES_spec_genes[, match(rownames(annotation_col), colnames(heatmap_counts_MES_spec_genes))]
heatmap_counts_ADRN_spec_genes <- heatmap_counts_ADRN_spec_genes[, match(rownames(annotation_col), colnames(heatmap_counts_ADRN_spec_genes))]


ensembl2symbol_annot <- ATAC_dds_results$results_all %>%
  dplyr::select(peak_id, gencode_gene_name)

#color.scheme <- rev(RColorBrewer::brewer.pal(8,"RdBu")) # generate the color scheme to use

color.scheme <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

ann_colors = list(
  cell_line = c( CLBM = "#005f73", NB10 = "#0a9396", NB6 = "#94d2bd", NB8 = "#DD3344", SH = "#FF9F1C"),
  cell_line_id = c(A = "#E9D8A6", M = "#D9BE6D")
)

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "Heatmap of signif. DAR promoters 3kb ADRN vs MES.",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

#ATAC_dds_results$results_signif[ATAC_dds_results$results_signif$log2FoldChange > 0,] 

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_ATAC.png"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_ATAC.eps"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")

# #now we can check how the peaks accessibility in MES and ADRN identified regions
# heatmap <- pheatmap::pheatmap(heatmap_counts_ADRN,
#                               main = "Heatmap of signif. DAR ADRN vs MES. promoters_3kb ADRN-specific",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN_promoters_3kb_enh_ATAC.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_ADRN_promoters_3kb_enh_ATAC.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# 
# heatmap <- pheatmap::pheatmap(heatmap_counts_MES,
#                               main = "Heatmap of signif. DAR ADRN vs MES. promoters 3kb MES-specific",
#                               scale = "row",
#                               annotation_col = annotation_col,
#                               annotation_colors = ann_colors,
#                               show_colnames = FALSE,
#                               show_rownames = FALSE,
#                               cluster_cols = FALSE,
#                               color = color.scheme,
#                               fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11
# 
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_promoters_3kb_enh_ATAC.png"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")
# ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_MES_promoters_3kb_enh_ATAC.eps"), 
#        plot=heatmap,
#        width = 18, height = 20, units = "cm")

heatmap <- pheatmap::pheatmap(heatmap_counts_MES_spec_genes,
                              main = "Heatmap of signif. DAR ADRN vs MES in promoters 3kb . MES-specific genes only (defined by RNA-seq)",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_MES-RNAseq-specific.png"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_MES-RNAseq-specific.eps"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")


heatmap <- pheatmap::pheatmap(heatmap_counts_ADRN_spec_genes,
                              main = "Heatmap of signif. DAR ADRN vs MES in promoters 3kb. ADRN-specific genes only (defined by RNA-seq)",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = FALSE,
                              cluster_cols = FALSE,
                              color = color.scheme,
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11

ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_ADRN-RNAseq-specific.png"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")
ggsave(filename = file.path(deg_dir, "heatmap_ADRN_vs_MES_promoters_3kb_ADRN-RNAseq-specific.eps"), 
       plot=heatmap,
       width = 18, height = 20, units = "cm")




# This part - we can identify binding sites of TFs and see if there they are differentially enriched between ADRN and MES lines
# Identifying TFBSs
# First, we load the PWM sets from JASPAR2022 DB.
# 9606 is a code for human
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["species"]] <- "9606"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- TFBSTools::getMatrixSet(JASPAR2022, opts)

ATAC_dds <- ATAC_dds_full
# Using different subsets of peaks to check 
subsets_to_run <- c("ALL", "ADRN_sh", "MES_sh", "g_ADRN", "g_MES")
subsets_to_run <- c("ALL")

for (sbst in subsets_to_run){
  sbst <-  c("ALL")
  # here we check if there are any peaks that have less than 5 reads across all samples.
  dim(ATAC_dds)
  counts_consensus_filt <- ATAC_dds_full[rowSums(assay(ATAC_dds_full)) > 5, ]
  dim(counts_consensus_filt)
  
  print(paste("procesing: ", sbst ))
  
  # if (sbst != "ALL"){
  #   if (sbst %in% c("ADRN_sh", "MES_sh")){
  #     counts_consensus_filt <- counts_consensus_filt[counts_consensus_filt@rowRanges@elementMetadata@listData[["sEnh_type"]] == sbst,]
  #   }else{
  #       counts_consensus_filt <- counts_consensus_filt[counts_consensus_filt@rowRanges@elementMetadata@listData[["gene_category"]] == sbst,]
  #   }
  # }

  
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
    dplyr::select(phenotype)
  pheatmap::pheatmap(as.dist(sample_cor), 
                     annotation_row = annotation_row_sampleCor, 
                     clustering_distance_rows = as.dist(1-sample_cor), 
                     clustering_distance_cols = as.dist(1-sample_cor) ) 
  
  #checkign the clusterrization using tSNE
  tsne_results <- chromVAR::deviationsTsne(chrom_access_deviations, threshold = 1.5, perplexity = 8, 
                                           what = "samples", shiny = FALSE)
  
  tsne_plots <- chromVAR::plotDeviationsTsne(chrom_access_deviations, tsne_results,
                                             sample_column = "phenotype", shiny = FALSE)
  
  diff_acc <- chromVAR::differentialDeviations(chrom_access_deviations, groups="phenotype", parametric = FALSE)
  motif_annot <- as.data.frame(SummarizedExperiment::rowData(chrom_access_deviations)) %>%
    tibble::rownames_to_column(var = "motif") %>%
    dplyr::select(motif, name)
  diff_acc_annot <- diff_acc %>%
    tibble::rownames_to_column(var = "motif") %>%
    dplyr::left_join(., motif_annot, by = "motif")
  head(diff_acc_annot) 
  
  diff_var <- chromVAR::differentialVariability(chrom_access_deviations, "phenotype", parametric = FALSE)
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
  ggsave(chrom_access_variability_plot, filename = file.path(deg_dir, paste0("ChromVAR_chrom_access_variability_corr_plot_", sbst, ".pdf")))
  
  BiocParallel::register(BiocParallel::SerialParam())
  
  chrom_access_variability_ord <- chrom_access_variability[order(chrom_access_variability$p_value), ]  # ordering based on p-value
  #chrom_access_variability_ord[1:10, ]
  
  #! to store
  chrom_access_variability_ord_results <- chrom_access_variability_ord %>%
    tibble::rownames_to_column(var="jaspar_id") %>%
    dplyr::left_join(., devZscores_df, by = "jaspar_id")
  
  #sum(chrom_access_variability_ord_results$p_value_adj < 0.05)
  
  #plotVariability(chrom_access_variability, use_plotly = FALSE) 
  
  message("Saving motif enrichment results")
  
  ntop <- 50
  topVariable <- chrom_access_variability_ord[1:ntop, ]
  topVariable <- tibble::rownames_to_column(topVariable, var="jaspar_id")
  
  devTop <- topVariable %>%
    dplyr::left_join(., devZscores_df, by = "jaspar_id")
  
  devTop <- devTop %>% mutate(long_name = paste0(jaspar_id, "_", name))
  
  devToPlot <- devTop %>%
    dplyr::select(long_name, CLBM_A_ATAC,  NB10_A_ATAC,  NB8_A_ATAC,  SH_A_ATAC, CLBM_M_ATAC, NB10_M_ATAC,  NB8_M_ATAC, SH_M_ATAC) %>%
    tibble::column_to_rownames(var = "long_name")
  #rownames(devToPlot) <- devTop[, 2]
  #library(pheatmap)
  
  annotCol_forHeatmap <- as.data.frame(colData(ATAC_dds_full)) %>%
    dplyr::select(phenotype)
  
  annotCol_forHeatmap_colors <- list(phenotype = c(A="#525252", M = "#fc4e2a"))
  
 # color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
  
  heatmap <- pheatmap::pheatmap(as.matrix(devToPlot),
                     #color = color.scheme,
                     color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
                     scale = "row",
                     #gaps_col = ph_col_annot%>% arrange(state)%$% table(state) %>% cumsum,
                     #gaps_row = ph_genes$clust %>% table() %>% cumsum,
                     cluster_cols = FALSE,
                     cluster_rows = TRUE,
                     show_colnames = TRUE,
                     show_rownames = TRUE,
                     annotation_col = annotCol_forHeatmap,
                     annotation_colors = annotCol_forHeatmap_colors,  
                     fontsize_row = 9,
                     fontsize = 9,
                     main = paste0("Heatmap of signif. Motifs ADRN vs MES. ", sbst, "- used"))
  
  ggsave(filename = file.path(deg_dir, paste0("heatmap_ADRN_vs_MES_Motifs", sbst, "- used", ".png")), 
         plot=heatmap,
         width = 18, height = 20, units = "cm")
  
  ggsave(filename = file.path(deg_dir, paste0("heatmap_ADRN_vs_MES_Motifs", sbst, "- used", ".eps")), 
         plot=heatmap,
         width = 18, height = 20, units = "cm")
}






## GSEA analysis #########
# Loading MsigDB geneset collections 
gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean=TRUE) 
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean=TRUE) 
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 

# Produce GSEA plots
MES_genes <- ATAC_dds_results$results_signif %>% dplyr::filter(log2FoldChange > 0) %>% pull(gencode_gene_name)
ADR_genes <- ATAC_dds_results$results_signif %>% dplyr::filter(log2FoldChange < 0) %>% pull(gencode_gene_name)


C5_GOBP_MES <- hypeR::hypeR(signature = MES_genes, 
                        genesets = gs_C5_GOBP, 
                        test = "hypergeometric", 
                        background = nrow(ATAC_dds))
C5_GOBP_ADR <- hypeR::hypeR(signature = ADR_genes, 
                        genesets = gs_C5_GOBP, 
                        test = "hypergeometric", 
                        background = nrow(ATAC_dds))
C5_GOBP_plot_MES <- hypeR::hyp_dots(C5_GOBP_MES, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: MES") +theme_bw()
C5_GOBP_plot_ADR <- hypeR::hyp_dots(C5_GOBP_ADR, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: ADR") +theme_bw()
C5_GOBP_plot_MES
C5_GOBP_plot_ADR

C2_kegg_MES <- hypeR::hypeR(signature = MES_genes, 
                        genesets = gs_C2_kegg, 
                        test="hypergeometric", 
                        background=nrow(ATAC_dds))
C2_kegg_ADR <- hypeR::hypeR(signature = ADR_genes, 
                        genesets = gs_C2_kegg, 
                        test="hypergeometric", 
                        background=nrow(ATAC_dds))
C2_kegg_plot_MES <- hypeR::hyp_dots(C2_kegg_MES, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="Kegg: MES") +theme_bw()
C2_kegg_plot_ADR <- hypeR::hyp_dots(C2_kegg_ADR, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="Kegg: ADR") +theme_bw()
C2_kegg_plot_MES
C2_kegg_plot_ADR

C2_reactome_MES <- hypeR::hypeR(signature = MES_genes, 
                        genesets = gs_C2_reactome, 
                        test="hypergeometric", 
                        background=nrow(ATAC_dds))
C2_reactome_ADR <- hypeR::hypeR(signature = ADR_genes, 
                        genesets = gs_C2_reactome, 
                        test="hypergeometric", 
                        background=nrow(ATAC_dds))
C2_reactome_plot_MES <- hypeR::hyp_dots(C2_reactome_MES, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="Reactome: MES") +theme_bw()
C2_reactome_plot_ADR <- hypeR::hyp_dots(C2_reactome_ADR, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="Reactome: ADR") +theme_bw()
C2_reactome_plot_MES

pdf(file = file.path(deg_dir, "MES_and_ADR_GSEA.pdf"), width = 20, height = 10)
C5_GOBP_plot_MES + C5_GOBP_plot_ADR
C2_kegg_plot_MES + C2_kegg_plot_ADR
C2_reactome_plot_MES + C2_reactome_plot_ADR
dev.off()


#optional saving to excel tables
#hypeR::hyp_to_excel(C5_GOBP, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOBP.xlsx"))
#hypeR::hyp_to_excel(C5_GOCC, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOCC.xlsx"))
#hypeR::hyp_to_excel(C5_GOMF, file_path=file.path(deg_dir, "MES_vs_ADRN_GSEA_C5_GOMF.xlsx"))





