library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

url <- "https://plus.figshare.com/ndownloader/files/51065297"
file_path <- "~/workspace/neuroblastoma/data/temp_data/Model.csv"

download.file(url = url, destfile = file_path)
Model_dataframe <- read.csv(file = file_path, header = TRUE, sep = ",")
Model_dataframe <- Model_dataframe %>% dplyr::select(ModelID, StrippedCellLineName, OncotreeSubtype)
Model_dataframe$MetaGroup <- ifelse(Model_dataframe$OncotreeSubtype == "Neuroblastoma", "Neuroblastoma", "Other")
Model_dataframe$MetaGroup <- ifelse(Model_dataframe$StrippedCellLineName %in% c("GIMEN", "NB69", "KPNSI9S"), 
                                    Model_dataframe$StrippedCellLineName,
                                    Model_dataframe$MetaGroup)

url <- "https://plus.figshare.com/ndownloader/files/51064667"
file_path <- "~/workspace/neuroblastoma/data/temp_data/ScreenGeneEffect.csv"
download.file(url = url, destfile = file_path)
Screen_Gene_Effect <- read.csv(file = file_path, header = TRUE, sep = ",")
colnames(Screen_Gene_Effect) <- str_replace_all(colnames(Screen_Gene_Effect), 
                                                          pattern = "\\.\\..*$", replacement = "")
Screen_Gene_Effect_compact <- Screen_Gene_Effect[, c("X", "WWTR1", "YAP1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "JUN", "FOS", "FOSL1", "FOSL2")]
colnames(Screen_Gene_Effect_compact)[1] <- "ModelID"
expression_long <- Screen_Gene_Effect_compact %>%
  pivot_longer(cols = -ModelID, names_to = "Gene", values_to = "GeneEffect") %>%
  left_join(Model_dataframe, by = "ModelID") %>%
  mutate(order_flag = ifelse(MetaGroup != "Other", 1, 0)) %>%
  arrange(order_flag) %>% 
  mutate(order_flag = ifelse(MetaGroup %in% c("GIMEN", "NB69", "KPNSI9S"), 1, 0)) %>%
  arrange(order_flag)

ggplot(expression_long, aes(x = Gene, y = GeneEffect, color = MetaGroup)) +
  geom_point() +
  scale_color_manual(
    values = c(
      "GIMEN" = "firebrick",
      "NB69" = "green",
      "KPNSI9S" = "yellow",
      "Neuroblastoma" = "navy",
      "Other" = "lightgrey"
    )
  ) +
  scale_alpha_manual(values = c(
      "GIMEN" = 1,
      "NB69" = 1,
      "KPNSI9S" = 1,
      "Neuroblastoma" = 1,
      "Other" = 0.2
  )) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Effect per Sample",
       x = "Sample (X)", y = "Gene Effect")


table((expression_long %>% filter(Gene == "FOS"))$MetaGroup)


######################################################################
######################################################################
######################################################################


# Visualisation of decon eda data (neuroblastoma plasticity clinical data)
library("ggplot2")
library("ggpubr")
library("GSVA")
library("pheatmap")
library("dplyr")
library("tidyr")
library(stringr)
library(org.Hs.eg.db)
# loading data ----
system("unzip ~/workspace/neuroblastoma/data/other/decon_eda_louis.RData.zip -d ~/workspace/neuroblastoma/data/other/")
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
                                       keys = c("CD44", "SOX10", "PLP1", "JUN", "VIM", "FOSL1", "FN1", "YAP1", "WWTR1", "FOSL2", "DBH", "PHOX2B", "PTPRC"),
                                       columns = c("ENSEMBL", "SYMBOL"),
                                       keytype = "SYMBOL")
genes_to_plot <- genes_to_plot %>% 
  filter(ENSEMBL != "ENSG00000262418")

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


### Visualize this as a boxplot
data_to_plot_subset_bxplt <- data_to_plot_subset %>%
  mutate(gene = row.names(.)) %>% 
  pivot_longer(cols = -gene,
               names_to = "sample",
               values_to = "expression"
  ) %>%
  mutate(
    group = case_when(
      str_starts(sample, "DTC_") ~ "DTC",
      str_starts(sample, "MNC_") ~ "MNC",
      str_starts(sample, "BMn_") ~ "BMn",
      TRUE ~ "Other")
  )
data_to_plot_subset_bxplt$group <- factor(data_to_plot_subset_bxplt$group, levels = c("DTC", "MNC", "BMn"))

data_to_plot_subset_bxplt %>% ggplot(aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ gene, scales = "free_y") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  theme_minimal()

data_to_plot_subset_bxplt %>% ggplot(aes(x = group, y = expression, fill = group)) +
  facet_wrap(~ gene, scales = "free_y") +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5,
               aes(ymin = ..y.., ymax = ..y..), color = "red", fatten = 0) +
  theme_minimal()



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

mes_adrn_gene_list$MES



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





####### check correlation to clinical parameters

RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
mes_adrn_gene_list <- list(
  MES = mes_adrn_gene_list[mes_adrn_gene_list$Term == "MES", "ensembl_id"],
  ADR = mes_adrn_gene_list[mes_adrn_gene_list$Term == "ADRN", "ensembl_id"]
)

load("~/workspace/neuroblastoma/data/other/decon_eda_louis.RData")
vst_counts_removedBatchEffect <- as.data.frame(vst_counts_removedBatchEffect) %>% 
  dplyr::select(starts_with("DTC"))

dx_rem <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect),
                               pattern = "(REL)|(DX)")
ssgsea_mes_adr_cellines <- GSVA::gsva(
  as.matrix(vst_counts_removedBatchEffect),
  mes_adrn_gene_list,
  method=c("ssgsea"),
  min.sz=1, 
  max.sz=Inf, 
  ssgsea.norm=TRUE, verbose=TRUE, parallel.sz=10)

selected_samples_names <- c("DTC_041_REL_2",
                            "DTC_061_DX_2",
                            "DTC_001_DX",
                            "DTC_010_REL",
                            "DTC_059_REL",
                            "DTC_010_DX",
                            "DTC_031_DX")

ssgsea_mes_adr_cellines <- as.data.frame(ssgsea_mes_adr_cellines) 
ssgsea_mes_adr_cellines <- as.data.frame(t(ssgsea_mes_adr_cellines))
ssgsea_mes_adr_cellines$MESARD_score <- ssgsea_mes_adr_cellines$MES - ssgsea_mes_adr_cellines$ADR
ssgsea_mes_adr_cellines$dx_rem <- dx_rem
ssgsea_mes_adr_cellines$sample_selector <- ifelse(row.names(ssgsea_mes_adr_cellines) %in% selected_samples_names,
                                                  "selected",
                                                  "therest")
ssgsea_mes_adr_cellines$YAP <- t(vst_counts_removedBatchEffect["ENSG00000137693",])
ssgsea_mes_adr_cellines$WWTR1 <- t(vst_counts_removedBatchEffect["ENSG00000018408",])
ssgsea_mes_adr_cellines$VIM <- t(vst_counts_removedBatchEffect["ENSG00000026025",])
ssgsea_mes_adr_cellines$JUN <- t(vst_counts_removedBatchEffect["ENSG00000177606",])



#### test - add more metadata
clinical_metadata <- read.csv2("~/workspace/neuroblastoma/data/other/Soeren_paper_Clinical_endpoints.csv")

indexes_to_use <- match(colnames(vst_counts_removedBatchEffect), clinical_metadata$OMICS_ID)
ssgsea_mes_adr_cellines <- cbind(ssgsea_mes_adr_cellines, 
                                 clinical_metadata[indexes_to_use, c("TP", "MNA", "Clinical_Relapse", "Dead_of_disease")]
)

ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(sample_selector, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression")

ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = sample_selector)) +
  geom_boxplot()


# checking expression between dx and rel
ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(dx_rem, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression") %>%
  filter(!is.na(dx_rem))

ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = dx_rem)) +
  geom_boxplot()

# checking expression between relapse and dx for Clinical_Relapse
ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(Clinical_Relapse, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression") %>%
  filter(!is.na(Clinical_Relapse))

ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = Clinical_Relapse)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Clinical_Relapse), method = "t.test") +
  labs(title = "Checking expression Clinical_Relapse param")

# checking expression Dead_of_disease param
ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(Dead_of_disease, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression") %>%
  filter(!is.na(Dead_of_disease))


ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = Dead_of_disease)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Dead_of_disease), method = "t.test") +
  labs(title = "Checking expression Dead_of_disease param")





rownames(ssgsea_mes_adr_cellines) == colnames(vst_counts_removedBatchEffect)

t.test(ssgsea_mes_adr_cellines$MESARD_score ~ ssgsea_mes_adr_cellines$sample_selector)
t.test(ssgsea_mes_adr_cellines$MESARD_score ~ ssgsea_mes_adr_cellines$dx_rem)
t.test(ssgsea_mes_adr_cellines$MESARD_score ~ ssgsea_mes_adr_cellines$Dead_of_disease)

t.test(ssgsea_mes_adr_cellines$YAP ~ ssgsea_mes_adr_cellines$dx_rem)
t.test(ssgsea_mes_adr_cellines$YAP ~ ssgsea_mes_adr_cellines$Dead_of_disease)

t.test(ssgsea_mes_adr_cellines %>% filter(dx_rem == "DX") %>% pull(MES),
       ssgsea_mes_adr_cellines %>% filter(dx_rem == "REL") %>% pull(MES))
t.test(ssgsea_mes_adr_cellines$MES ~ ssgsea_mes_adr_cellines$dx_rem)

ssgsea_mes_adr_cellines$dx_rem <- as.factor(ssgsea_mes_adr_cellines$dx_rem) 
model <- glm(dx_rem ~ YAP, data = ssgsea_mes_adr_cellines, family = "binomial")
summary(model)

### Deal with repeating patients - take the average
ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines %>% 
  mutate(Sample_ID = sub("(_REL.*|_DX.*)", "", row.names(.)))  # Extract everything before the first underscore
ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines %>%
  group_by(Sample_ID) %>%
  summarize(
    across(where(is.numeric), mean, na.rm = TRUE),  # Average numeric columns
    across(where(is.character), dplyr::first),            # Retain first non-numeric value
    across(where(is.factor), dplyr::first),            # Retain first non-numeric value
    .groups = "drop"
  )

# removing "REL" samples, Sabine told that it doesn't make sense to include them here
ssgsea_mes_adr_cellines <- ssgsea_mes_adr_cellines %>% 
  filter(dx_rem == "DX")

# checking expression between relapse and dx for Clinical_Relapse
ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(Clinical_Relapse, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression") %>%
  filter(!is.na(Clinical_Relapse))

ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = Clinical_Relapse)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Clinical_Relapse), method = "t.test") +
  labs(title = "Checking expression Clinical_Relapse param REPEATING Samples removed")


# checking expression Dead_of_disease param
ssgsea_mes_adr_cellines_bxplot <- ssgsea_mes_adr_cellines %>% dplyr::select(Dead_of_disease, YAP, WWTR1, JUN, VIM)
ssgsea_mes_adr_cellines_bxplot <- tidyr::pivot_longer(ssgsea_mes_adr_cellines_bxplot, 
                                                      cols = c(YAP, WWTR1, JUN, VIM), 
                                                      names_to = "gene", 
                                                      values_to = "Expression") %>%
  filter(!is.na(Dead_of_disease))

ggplot(ssgsea_mes_adr_cellines_bxplot, aes(x = gene, y = Expression, fill = Dead_of_disease)) +
  geom_boxplot() +
  stat_compare_means(aes(group = Dead_of_disease), method = "t.test") +
  labs(title = "Checking expression Dead_of_disease param REPEATING Samples removed")



# Split the data into two groups - high and low YAP TAZ expression (1st and 4th quartiles)
# And check if the expression is changing between the dx_rem is different
load("~/workspace/neuroblastoma/data/other/decon_eda_louis.RData")
vst_counts_removedBatchEffect <- as.data.frame(vst_counts_removedBatchEffect) %>% 
  dplyr::select(starts_with("DTC"))

dx_rem <- stringr::str_extract(string = colnames(vst_counts_removedBatchEffect),
                               pattern = "(REL)|(DX)")

#YAP1 ENSG00000137693
#WWTR1 ENSG00000018408

quantiles <- quantile(t(vst_counts_removedBatchEffect)[,"ENSG00000137693"], 
                      probs = c(0.25, 0.75))

chi_square_data_table <- data.frame(DX_rem = dx_rem)

chi_square_data_table$YAP_expression <- t(vst_counts_removedBatchEffect["ENSG00000137693", ])


chi_square_data_table <- chi_square_data_table %>%
  mutate(
    YAP_group = case_when(
      YAP_expression <= quantiles[1] ~ "Low",
      YAP_expression >= quantiles[2] ~ "High",
      TRUE ~ "Middle"
    )
  )

chi_square_data_table <- chi_square_data_table %>% 
  filter(!is.na(DX_rem)) %>%
  filter(YAP_group %in% c("High", "Low"))

contingency_table <- table(chi_square_data_table$YAP_group, chi_square_data_table$DX_rem)
# View the table
print(contingency_table)

chisq_test <- chisq.test(contingency_table)
print(chisq_test)

fisher_test <- fisher.test(contingency_table)
print(fisher_test)








######################################################################
######################################################################
######################################################################

# Drug scoring 
# Code produced by Sören
library(tidyverse)

data_dir <- "~/workspace/neuroblastoma/data/other/"
results_dir <- "~/workspace/neuroblastoma/results/other/"

# Function to calculate AUC scores for drug screening data
calculate_auc_scores <- function(viabilities_file, output_file) {
  # Read and prepare data
  viabilities <- read.csv(file.path(data_dir, viabilities_file))
  viabilities <- viabilities %>% dplyr::rename("ADR" = NumberOfCells_ADR,
                                        "MES" = NumberOfCells_MES)
  controls <- viabilities %>% filter(Drug == "DMSO")
  
  ref_ADR <- mean(controls$ADR)
  ref_MES <- mean(controls$MES)
  
  viabilities <- viabilities %>% mutate(ADR_rel = ADR/ref_ADR, MES_rel = MES/ref_MES)
  
  druglist <- unique(viabilities$Drug)
  druglist <- druglist[druglist != "DMSO"]
  
  scores <- data.frame(drug = "", score_ADR = NA, score_MES = NA)
  
  for (drug in druglist) {
    temp2 <- dplyr::filter(viabilities, Drug == drug)
    concentrations <- sort(unique(temp2$Concentration))
    
    adr_means <- c()
    mes_means <- c()
    for (i in concentrations) {
      # Calculate mean values for viabilities at concentrations
      adr_means <- append(adr_means, c(mean(filter(temp2, Concentration == i)$ADR_rel)))
      mes_means <- append(mes_means, c(mean(filter(temp2, Concentration == i)$MES_rel)))
    }
    
    # Calculate AUC using trapezoidal rule
    auc_adr <- 0
    auc_mes <- 0
    l <- length(concentrations)
    for (i in c(2:(l))) {
      auc_adr <- auc_adr + ((1/(l-1)) * 0.5 * (adr_means[i] + adr_means[i-1]))
      auc_mes <- auc_mes + ((1/(l-1)) * 0.5 * (mes_means[i] + mes_means[i-1]))
    }
    
    # Convert to area above curve (drug efficacy score)
    score_adr <- (1-auc_adr)
    score_mes <- (1-auc_mes)
    
    # Add to scores dataframe
    scores <- add_row(scores, drug = drug, score_ADR = score_adr, score_MES = score_mes)
  }
  
  # Write results
  write.csv(scores, file = file.path(results_dir, output_file))
  return(scores)
}

# Calculate AUC scores for all cell lines
SK_N_SH_scores <- calculate_auc_scores("SK-N-SH_viabilities.csv", "SK-N-SH_reanalysis_v3.csv")
STA_NB_8_scores <- calculate_auc_scores("STA-NB-8_viabilities.csv", "STA-NB-8_reanalysis_v3.csv")
STA_NB_10_scores <- calculate_auc_scores("STA-NB-10_viabilities.csv", "STA-NB-10_reanalysis_v3.csv")
CLB_Ma_scores <- calculate_auc_scores("CLB-Ma_viabilities.csv", "CLB-Ma_reanalysis_v3.csv")

# Read all scores lists
SST055 <- read.csv(file.path(results_dir, "SK-N-SH_reanalysis_v3.csv"))
SST062 <- read.csv(file.path(results_dir, "STA-NB-10_reanalysis_v3.csv"))
SST064 <- read.csv(file.path(results_dir, "STA-NB-8_reanalysis_v3.csv"))
SST070 <- read.csv(file.path(results_dir, "CLB-Ma_reanalysis_v3.csv"))

#List of experiments with names
experiment_names <- c("SST055", "SST062", "SST064", "SST070")
experiments <- list(SST055, SST062, SST064, SST070)
names(experiments) <- experiment_names

#loop through all experiments for positive control correction
for (i in seq_along(experiments)) {
  # Access the dataframe
  df <- experiments[[i]]
  # Perform the calculations and update the dataframe
  df <- df %>%
    arrange(drug) %>%
    mutate(adj_score_ADR = if_else(score_ADR < 0, 0, score_ADR)) %>%
    mutate(adj_score_MES = if_else(score_MES < 0, 0, score_MES)) %>%
    mutate(adj_anti_MES = adj_score_MES - adj_score_ADR)
  
  # Save the modified dataframe back to the original variable
  assign(experiment_names[i], df)
}

write.csv(SST055, file = file.path(results_dir,"SK-N-SH_scores.csv"))
write.csv(SST062, file = file.path(results_dir,"STA-NB-10_scores.csv"))
write.csv(SST064, file = file.path(results_dir,"STA-NB-8_scores.csv"))
write.csv(SST070, file = file.path(results_dir,"CLB-Ma_scores.csv"))

#remove DMSO and empty row, outlier obatoclax (only results from 2 datapoints),
# and removal of retinoids: isotretinoin, tretinoin, bexarotene, alitretinoin (affects VIM staining), 
# and Verteporfin (wrong concentration printed)
# Define excluded drugs once
excluded_drugs <- c("DMSO", "", "Obatoclax", "Isotretinoin", "Tretinoin", "Alitretinoin", "Bexarotene", "Verteporfin")
# Apply filter to all datasets
SST055 <- SST055 %>% filter(!drug %in% excluded_drugs)
SST062 <- SST062 %>% filter(!drug %in% excluded_drugs)
SST064 <- SST064 %>% filter(!drug %in% excluded_drugs)
SST070 <- SST070 %>% filter(!drug %in% excluded_drugs)

#SK-N-SH top 10 top/bottom
top_hits <- top_n(SST055, 10, adj_anti_MES)
bot_hits <- top_n(SST055, -10, adj_anti_MES)
tb_hits <- rbind(top_hits, bot_hits)
ggplot(tb_hits, aes(x = adj_anti_MES, y = reorder(drug, adj_anti_MES), fill = adj_anti_MES)) +
  geom_bar(stat = "summary", fun = "mean") + 
  scale_fill_gradient2(name = "anti-MES", 
                       high = "#EB363C", 
                       low = "blue", mid = "grey95") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8, color = "black", face = "bold"), aspect.ratio = 1/1) + 
  labs(x = "anti-MES score", y = "")
ggsave(file.path(results_dir,"SK-N-SH_individual_drugs_tb10_v3.pdf"))

#STA-NB-10 top hits
top_hits <- SST062 %>% filter(adj_anti_MES > 0.00)
ggplot(top_hits, aes(x = adj_anti_MES, 
                     y = reorder(drug, adj_anti_MES), 
                     fill=adj_anti_MES)) + 
  geom_bar(stat = "summary", fun = "mean") + 
  scale_fill_gradient2(name = "anti-MES",
                       high = "#EB363C", 
                       low="blue",
                       mid="grey95") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 8, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave(file.path(results_dir, "NB-10_individual_drugs_tophits_v3.pdf"))

#STA-NB-8 top hits
top_hits <- SST064 %>% filter(adj_anti_MES > 0.19)
ggplot(top_hits,
       aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + 
  geom_bar(stat = "summary", fun = "mean") + 
  scale_fill_gradient2(name = "anti-MES", 
                       high="#EB363C",
                       low="blue",
                       mid="grey95") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 6, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave(file.path(results_dir, "NB-8_individual_drugs_tophits_v3.pdf"))

#CLB-Ma top hits
top_hits <- SST070 %>% filter(adj_anti_MES > 0.19)
ggplot(top_hits,
       aes(x=adj_anti_MES, y= reorder(drug, adj_anti_MES), fill=adj_anti_MES)) + 
  geom_bar(stat = "summary", fun = "mean") + 
  scale_fill_gradient2(name = "anti-MES",
                       high = "#EB363C",
                       low = "blue",
                       mid = "grey95") +
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 3.3, color = "black", face = "bold")) + 
  labs(x = "anti-MES score", y = "")
ggsave(file.path(results_dir, "CLB-Ma_individual_drugs_tophits_v3.pdf"))


# combine anti-MES values in one data frame and save
all_lines_anti_MES <- data.frame(drug = SST055$drug, 
                                 SK_N_SH = SST055$adj_anti_MES, 
                                 CLB_Ma = SST070$adj_anti_MES, 
                                 STA_NB_8 = SST064$adj_anti_MES, 
                                 STA_NB_10 = SST062$adj_anti_MES)
write.csv(all_lines_anti_MES, file.path(results_dir, "all_lines_anti-MES_v3.csv"))

# combine AUC_ADR values in one data frame and save
all_lines_AUC_ADR <- data.frame(drug = SST055$drug, 
                                SK_N_SH = SST055$adj_score_ADR, 
                                CLB_Ma = SST070$adj_score_ADR, 
                                STA_NB_8 = SST064$adj_score_ADR, 
                                STA_NB_10 = SST062$adj_score_ADR)
write.csv(all_lines_AUC_ADR, file.path(results_dir, "all_lines_AUC_ADR_v3.csv"))

# combine AUC_MES values in one data frame and save
all_lines_AUC_MES <- data.frame(drug = SST055$drug, 
                                SK_N_SH = SST055$adj_score_MES, 
                                CLB_Ma = SST070$adj_score_MES, 
                                STA_NB_8 = SST064$adj_score_MES, 
                                STA_NB_10 = SST062$adj_score_MES)
write.csv(all_lines_AUC_MES, file.path(results_dir, "all_lines_AUC_MES_v3.csv"))

#combined plot anti-MES all cell lines
all_lines <- read.csv(file.path(results_dir, "all_lines_anti-MES_v3.csv"))

numeric_cols <- sapply (all_lines, is.numeric)
numeric_cols["X"] <- FALSE  # Exclude column "X" from the condition
to_keep <- apply(all_lines[, numeric_cols], 1, function(row){
  all(row > 0, na.rm = TRUE) || all(row < 0, na.rm = TRUE)
})

only_consistent <- all_lines[to_keep,]

individual_long <- only_consistent %>%
  select(drug, SK_N_SH, CLB_Ma, STA_NB_8, STA_NB_10) %>%
  gather(key = "cell_line", value = "value", -drug)

summary_stats <- individual_long %>%
  group_by(drug) %>%
  summarise(
    avg_anti_MES = mean(value),
    sem = sd(value)/sqrt(4),
    .groups = "drop"
  )

top_hits <- summary_stats %>% top_n(-10, avg_anti_MES)
bot_hits <- summary_stats %>% top_n(10, avg_anti_MES)
tb_summary <- bind_rows(top_hits, bot_hits)

tb_drugs <- tb_summary$drug
individual_long_tb <- individual_long %>%
  filter(drug %in% tb_drugs)

tb_summary$drug <- factor(tb_summary$drug, levels = unique(tb_summary$drug))
individual_long_tb$drug <- factor(individual_long_tb$drug, 
                                  levels = levels(tb_summary$drug))

tb_summary$drug <- factor(tb_summary$drug, levels = tb_summary %>% arrange(avg_anti_MES) %>% pull(drug))
individual_long_tb$drug <- factor(individual_long_tb$drug, 
                                  levels = levels(tb_summary$drug))

# Definition custom shapes for geom_jitter
shape_values <- c(
  "SK_N_SH" = 16,    
  "CLB_Ma" = 17,     
  "STA_NB_8" = 15,   
  "STA_NB_10" = 18   
)


ggplot() +
  geom_bar(data = tb_summary, aes(y = drug, x = avg_anti_MES, fill = avg_anti_MES), stat = "identity") +
  geom_errorbarh(data = tb_summary, 
                 aes(y = drug, xmin = avg_anti_MES - sem, xmax = avg_anti_MES + sem), 
                 height = 0.3, 
                 linewidth = 0.2, 
                 color = "black") +
  scale_shape_manual(values = shape_values) +
  geom_jitter(data = individual_long_tb, 
              aes(y = drug, x = value, shape = cell_line), 
              width = 0.01, 
              height = 0.01, 
              color = "black", 
              size = 0.75) +
  scale_fill_gradient2(name = "anti_MES", high = "#EB363C", low = "blue", mid = "grey95") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 8, color = "black", face = "bold"),
        aspect.ratio = 1 / 1) +
  labs(y = "Drug", x = "anti-MES score")

ggsave(file.path(results_dir, "combined_drugs_tb_10_datapoints_symbols_consistent_error_bars_v3.pdf"))
