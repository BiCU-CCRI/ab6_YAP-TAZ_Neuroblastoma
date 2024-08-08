library("GSVA")
library("pheatmap")
library("dplyr")
library(stringr)
library(org.Hs.eg.db)
# loading data ----
system("unzip ~/workspace/neuroblastoma/data/other/decon_eda_louis.RData.zip")
load("~/workspace/neuroblastoma/data/other/decon_eda_louis.RData")

# row_annotation_table <- data.frame(Symbol = genes_to_plot$SYMBOL, row.names = genes_to_plot$ENSEMBL)
# data_to_plot <- as.data.frame(vst_counts_removedBatchEffect) %>% filter(row.names(.) %in% genes_to_plot$ENSEMBL)
# row.names(data_to_plot) <- genes_to_plot$SYMBOL[match(row.names(data_to_plot), genes_to_plot$ENSEMBL)]
# data_to_plot %>% 
#   dplyr::select(starts_with(c("DTC", "MNC", "BMn"))) 
# %>%
#   dplyr::select(-matches("MNC_0[1-9].")) %>%
#   dplyr::select(-matches("BMn_[6-9]."))







# # load genes from RNA-seq MES/ADR-specific genes
# RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/20240515/cell_type_MES_vs_ADR.xlsx", sheet = 1)
# mes_adrn_gene_list <- RNA_SEQ_data %>%
#   mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
# mes_adrn_gene_list <- list(
#   MES = mes_adrn_gene_list[mes_adrn_gene_list$Term == "MES", "ensembl_id"],
#   ADR = mes_adrn_gene_list[mes_adrn_gene_list$Term == "ADRN", "ensembl_id"]
# )
# 
# signatures_raw <- readr::read_tsv(file = "~/workspace/neuroblastoma/data_soren/legacy_data_mes_adr_ncc_noradr_genesets.tsv", skip_empty_rows=TRUE) # do not add NA for empty
# signatures_list <- list(Mesenchymal_groen  = signatures_raw$Mesenchymal[!is.na(signatures_raw$Mesenchymal)],
#                         Adrenergic_groen = signatures_raw$Adrenergic[!is.na(signatures_raw$Adrenergic)])
# signatures_list$Mesenchymal_groen <- AnnotationDbi::select(hs,
#                                                            keys = signatures_list$Mesenchymal_groen,
#                                                            columns = c("ENSEMBL", "SYMBOL"),
#                                                            keytype = "SYMBOL") %>% 
#   filter(!is.na(ENSEMBL)) %>%
#   pull(ENSEMBL)
# signatures_list$Adrenergic_groen <- AnnotationDbi::select(hs,
#                                                            keys = signatures_list$Adrenergic_groen,
#                                                            columns = c("ENSEMBL", "SYMBOL"),
#                                                            keytype = "SYMBOL") %>% 
#   filter(!is.na(ENSEMBL)) %>%
#   pull(ENSEMBL)
# signatures_list <- c(signatures_list, mes_adrn_gene_list)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# the real plot
hs <- org.Hs.eg.db
genes_to_plot <- AnnotationDbi::select(hs,
                                       keys = c("CD44", "SOX10", "PLP1", "JUN", "VIM", "FOSL1", "FN1", "YAP1", "WWTR1", "FOSL2", "DBH", "PHOX2B"),
                                       columns = c("ENSEMBL", "SYMBOL"),
                                       keytype = "SYMBOL")

load("~/workspace/neuroblastoma/data/other/decon_eda_louis.RData")
vst_counts_removedBatchEffect_subset <- as.data.frame(vst_counts_removedBatchEffect) %>% 
  dplyr::select(starts_with(c("DTC", "MNC", "BMn")))

tmp_ids_subset <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect_subset), 
                                pattern = "(DTC)|(MNC)|(BMn)")
dx_rem <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect_subset),
                               pattern = "(REL)|(DX)")
heatmap_col_annot_DTC_subset <- data.frame(id = tmp_ids_subset, dx_rem = dx_rem)
row.names(heatmap_col_annot_DTC_subset) <- colnames(vst_counts_removedBatchEffect_subset)

row_annotation_table <- data.frame(Symbol = genes_to_plot$SYMBOL, 
                                   row.names = genes_to_plot$ENSEMBL)
data_to_plot_subset <- vst_counts_removedBatchEffect_subset %>% 
  filter(row.names(.) %in% genes_to_plot$ENSEMBL)
row.names(data_to_plot_subset) <- genes_to_plot$SYMBOL[match(row.names(data_to_plot_subset), genes_to_plot$ENSEMBL)]

DTC_enrichm_score <- read.csv2("~/workspace/neuroblastoma/data/other/DTC_enrichment_efficiency.csv",
                               header = TRUE)
heatmap_col_annot_DTC_subset <- merge.data.frame(x = DTC_enrichm_score, 
                                                 y = heatmap_col_annot_DTC_subset, 
                                                 by.y = "row.names", 
                                                 by.x = "OMICS_ID", 
                                                 all.y = TRUE)                      
row.names(heatmap_col_annot_DTC_subset) <- heatmap_col_annot_DTC_subset$OMICS_ID
heatmap_col_annot_DTC_subset <- heatmap_col_annot_DTC_subset %>% 
  dplyr::select(-OMICS_ID)
heatmap_col_annot_DTC_subset[is.na(heatmap_col_annot_DTC_subset$percent_DTC_AE), "percent_DTC_AE"] <- 0
heatmap_col_annot_DTC_subset[heatmap_col_annot_DTC_subset$id != "DTC", "dx_rem"] <- "NA"
heatmap_col_annot_DTC_subset[is.na(heatmap_col_annot_DTC_subset$dx_rem), "dx_rem"] <- "NA"

colorRampPalette(c("white", "green"))(26)
color_vector <- c(colorRampPalette(c("darkgrey"))(1), 
                  colorRampPalette(c("white", "forestgreen"))(26))
names(color_vector) <- c("NAA", sort(unique(heatmap_col_annot_DTC_subset$percent_DTC_AE)))
dx_rem_col_vector <- c("grey", "chocolate1", "deepskyblue1")
names(dx_rem_col_vector) <- c("NA", "DX", "REL")
ann_colors = list(
  percent_DTC_AE = c(color_vector),
  dx_rem = c(dx_rem_col_vector)
)

# data_to_plot_subset <- data_to_plot_subset %>%
#   dplyr::select(-matches("MNC_0[1-9].")) %>%
#   dplyr::select(-matches("BMn_[6-9]."))

pheatmap(data_to_plot_subset,
         scale = "row",
         annotation_col = heatmap_col_annot_DTC_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_colors = ann_colors,
         show_colnames = TRUE
)


# now - the GSVA data
# load genes from RNA-seq MES/ADR-specific genes
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
mes_adrn_gene_list <- list(
  MES = mes_adrn_gene_list[mes_adrn_gene_list$Term == "MES", "ensembl_id"],
  ADR = mes_adrn_gene_list[mes_adrn_gene_list$Term == "ADRN", "ensembl_id"]
)

signatures_raw <- readr::read_tsv(file = "~/workspace/neuroblastoma/data/other/legacy_data_mes_adr_ncc_noradr_genesets.tsv", skip_empty_rows=TRUE) # do not add NA for empty
signatures_list <- list(Mesenchymal_groen  = signatures_raw$Mesenchymal[!is.na(signatures_raw$Mesenchymal)],
                        Adrenergic_groen = signatures_raw$Adrenergic[!is.na(signatures_raw$Adrenergic)])
signatures_list$Mesenchymal_groen <- AnnotationDbi::select(hs,
                                                           keys = signatures_list$Mesenchymal_groen,
                                                           columns = c("ENSEMBL", "SYMBOL"),
                                                           keytype = "SYMBOL") %>% 
  filter(!is.na(ENSEMBL)) %>%
  pull(ENSEMBL)
signatures_list$Adrenergic_groen <- AnnotationDbi::select(hs,
                                                          keys = signatures_list$Adrenergic_groen,
                                                          columns = c("ENSEMBL", "SYMBOL"),
                                                          keytype = "SYMBOL") %>% 
  filter(!is.na(ENSEMBL)) %>%
  pull(ENSEMBL)
signatures_list <- c(signatures_list, mes_adrn_gene_list)





vst_counts_removedBatchEffect <- as.data.frame(vst_counts_removedBatchEffect) %>% 
  dplyr::select(starts_with("DTC") | starts_with("MNC") | starts_with("BMn"))

tmp_ids <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect_subset), 
                                pattern = "(DTC)|(MNC)|(BMn)")
dx_rem <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect_subset),
                               pattern = "(REL)|(DX)")

heatmap_col_annot_DTC <- data.frame(id = tmp_ids, dx_rem = dx_rem)
heatmap_col_annot_DTC[heatmap_col_annot_DTC$id != "DTC", "dx_rem"] <- NA

row.names(heatmap_col_annot_DTC) <- colnames(vst_counts_removedBatchEffect)


ssgsea_mes_adr_cellines <- GSVA::gsva(
  as.matrix(vst_counts_removedBatchEffect),
                                      signatures_list,
                                      method=c("ssgsea"),
                                      min.sz=1, 
                                      max.sz=Inf, 
                                      ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)

ssgsea_mes_adr_cellines <- data.frame(ssgsea_mes_adr_cellines) 

# FACS enrichment 
# 
# DTC_enrichm_score <- read.csv2("~/workspace/neuroblastoma/resources/DTC_enrichment_efficiency.csv",
#                                header = TRUE)
# heatmap_col_annot_DTC_subset <- data.frame(id = tmp_ids_subset)
# row.names(heatmap_col_annot_DTC_subset) <- colnames(vst_counts_removedBatchEffect_subset)
# heatmap_col_annot_DTC_subset <- merge.data.frame(x = DTC_enrichm_score, y = heatmap_col_annot_DTC_subset, by.y = "row.names", by.x = "OMICS_ID", all.y = TRUE)                      
# row.names(heatmap_col_annot_DTC_subset) <- heatmap_col_annot_DTC_subset$OMICS_ID
# heatmap_col_annot_DTC_subset <- heatmap_col_annot_DTC_subset %>% dplyr::select(-OMICS_ID)
# 
# colorRampPalette(c("white", "green"))(26)
# color_vector <- c(colorRampPalette(c("darkgrey"))(1), colorRampPalette(c("white", "forestgreen"))(26))
# names(color_vector) <- c("NAA", sort(unique(heatmap_col_annot_DTC_subset$percent_DTC_AE)))
# heatmap_col_annot_DTC_subset[is.na(heatmap_col_annot_DTC_subset$percent_DTC_AE), "percent_DTC_AE"] <- "NAA"
# ann_colors = list(
#   percent_DTC_AE = c(color_vector)
# )

pheatmap(ssgsea_mes_adr_cellines,
         scale = "row",
         annotation_col = heatmap_col_annot_DTC_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_colors = ann_colors,
         show_colnames = TRUE
)
