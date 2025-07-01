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


### Visualize this as boxplot
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





####### check correlation to clinical parameter

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




