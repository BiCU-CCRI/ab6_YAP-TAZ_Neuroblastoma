######
## Integration of ATAC-seq data with public HiC data.
######

#### Loading libraries ####
import::from(stringr, str_extract, str_detect)
library(dplyr)
import::from(.from = magrittr, "%>%") 
import::from(IRanges, IRanges)
import::from(GenomicRanges,GRanges)
library(DESeq2)
import::from(rtracklayer, liftOver, import.chain)
import::from(openxlsx, addWorksheet, writeData, saveWorkbook)
import::from(openxlsx2, read_xlsx)
import::from(.from = JASPAR2022, JASPAR2022)
import::from(.from = BSgenome.Hsapiens.UCSC.hg38, BSgenome.Hsapiens.UCSC.hg38)
import::from(.from = TFBSTools, getMatrixSet)
library(BSgenome.Hsapiens.UCSC.hg38)
import::from(chromVAR, addGCBias)
import::from(.from = RColorBrewer, brewer.pal)
import::from(.from = "~/workspace/neuroblastoma/resources/utilityScripts.R", "generatePCA", "extract_results_DDS", "meanExprsPerGroup", "extract_results_DDS_HIC")



##################
#filtering potential enhancer elements of YAP1

scale2_bed <- read.table("~/workspace/neuroblastoma/data_public/Hi-C/processed/scale2_bed.csv", sep = ";", header = T )
scale2_bed <- scale2_bed%>%select(-one_of("X"))
scale2_bed <- GRanges(scale2_bed)

enh_peak_ids <- subsetByOverlaps(ATAC_dds@rowRanges, scale2_bed,ignore.strand = TRUE)

ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)
vsd <- vst(ATAC_dds, blind = F)

ATAC_dds <- ATAC_dds[rowData(ATAC_dds)$peak_id %in% enh_peak_ids$peak_id,]
vsd <- vsd[rowData(vsd)$peak_id %in% enh_peak_ids$peak_id,]

resultsNames(ATAC_dds)
coeff_name <- 'phenotype_M_vs_A'
cond_numerator <-  "A"
cond_denominator <-  "M"
cond_variable <-  "phenotype"
padj_cutoff = 0.05
log2FC_cutoff = 0.58

ATAC_dds_results <- extract_results_DDS(dds_object = ATAC_dds,
                                        coeff_name = coeff_name,
                                        cond_numerator = cond_numerator,
                                        cond_denominator = cond_denominator,
                                        cond_variable = cond_variable,
                                        padj_cutoff = 0.05,
                                        log2FC_cutoff = 0.58)

metadata_heatmap <- as.data.frame(colData(ATAC_dds))
heatmap_counts <- SummarizedExperiment::assay(vsd)
rownames(heatmap_counts) <- vsd@rowRanges$PeakId

heatmap_counts <- SummarizedExperiment::assay(vsd) 
rownames(heatmap_counts) <- vsd@rowRanges$peak_id
# removed xperiment batch effects! ? use experimen+cell_line removed???

heatmap_counts <- heatmap_counts[rownames(heatmap_counts) %in% ATAC_dds_results$results_signif$peak_id,]

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, phenotype) %>% 
  dplyr::arrange(phenotype, cell_line)

heatmap_counts<- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]

ensembl2symbol_annot <- ATAC_dds_results$results_all %>%
  dplyr::select(peak_id, gencode_gene_name)

#color.scheme <- rev(RColorBrewer::brewer.pal(8,"RdBu")) # generate the color scheme to use

#Alternative color schemes
#color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012", "#9B2226")
color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
#color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")

ann_colors = list(
  cell_line = c( CLBM = "#005f73", NB10 = "#0a9396", NB8 = "#DD3344", SH = "#FF9F1C"),
  cell_line_id = c(A = "#E9D8A6", M = "#D9BE6D")
)

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "Heatmap of signif. DAR ADRN vs MES.",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = T,
                              cluster_cols = FALSE,
                              #cluster_rows = counts_deg_ord_row_cor_hclust,
                              color = color.scheme,
                              #cutree_rows = 5,
                              #gaps_col = gap_cols,
                              #labels_row = make_bold_names(vsd_tp2_vs_tp1_inclEXP9_7_counts_deg_ord_symbols, rownames, rc_names = names(genes_of_interest)),
                              fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11













#
##############
# Analysis of the peaks that overlaps with enhancers that were identified in:
#  https://sci-hub.ru/https://www.nature.com/articles/ncb3216
#

MDA_enhancers <- read.table("~/workspace/neuroblastoma/resources/MDA-MB-231_enhancers.csv", header=TRUE, sep = ";")
MDA_enhancers_IRanges <- GRanges(MDA_enhancers)

# run liftOver - it's necessary for converting hg19 to hg38 coordinates
chainObject <- import.chain("./neuroblastoma/resources/hg19ToHg38.over.chain")
MDA_enhancers_IRanges <- GRanges(as.data.frame(liftOver(MDA_enhancers_IRanges, chainObject)))

ATAC_dds <- readRDS("~/workspace/neuroblastoma/RDSs/ATAC_dds.RDS")

subset_atac <- subsetByOverlaps(rowRanges(ATAC_dds), MDA_enhancers_IRanges)
MDA_ranges <- subsetByOverlaps(MDA_enhancers_IRanges,rowRanges(ATAC_dds))
#subset_atac$target_gene <- subsetByOverlaps(MDA_enhancers_IRanges,rowRanges(ATAC_dds))$TARGET.GENE

result_atac <- GRanges()
tmp_atac <- GRanges()
for (range_number in 1:length(MDA_ranges)){
  tmp_atac <- subsetByOverlaps(subset_atac, MDA_ranges[range_number])
  tmp_atac$ENH <- MDA_ranges[range_number]$TARGET.GENE
  result_atac <- append(result_atac, tmp_atac)
  print(result_atac)
}



result_atac_bckup <- result_atac
table(result_atac$peak_id)
#remove duplicates
result_atac <- result_atac[!duplicated(result_atac$peak_id),]
ATAC_dds <- ATAC_dds[which(rowData(ATAC_dds)$peak_id %in% result_atac$peak_id),]
result_atac_reduced <- as.data.frame(result_atac) %>% select(ENH, peak_id)
result_atac <- merge(as.data.frame(rowData(ATAC_dds)), result_atac_reduced, by = ("peak_id"), suffixes = c("", ""))

#test that the order is correct. if 0 then it's good
sum(rowData(ATAC_dds)$peak_id != result_atac$peak_id)

rowData(ATAC_dds) <- result_atac

# optional - remove all entries that are not associated with a particular gene
#ATAC_dds <- ATAC_dds[rowData(ATAC_dds)$ENH != "",]

ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)

vsd <- vst(ATAC_dds, blind = F)

resultsNames(ATAC_dds)
coeff_name <- 'phenotype_M_vs_A'
cond_numerator <-  "A"
cond_denominator <-  "M"
cond_variable <-  "phenotype"
padj_cutoff = 0.05
log2FC_cutoff = 0.58

ATAC_dds_results <- extract_results_DDS_HIC(dds_object = ATAC_dds,
                                            coeff_name = coeff_name,
                                            cond_numerator = cond_numerator,
                                            cond_denominator = cond_denominator,
                                            cond_variable = cond_variable,
                                            padj_cutoff = 0.05,
                                            log2FC_cutoff = 0.58)

XLSX_OUT <- openxlsx::createWorkbook()
openxlsx::addWorksheet(XLSX_OUT, "results_signif")
openxlsx::addWorksheet(XLSX_OUT, "de_details")
openxlsx::addWorksheet(XLSX_OUT, "results_all")

openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_signif, sheet = "results_signif")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$de_details, sheet = "de_details")
openxlsx::writeData(XLSX_OUT, x = ATAC_dds_results$results_all, sheet = "results_all")

openxlsx::saveWorkbook(XLSX_OUT, file.path("~/workspace/neuroblastoma/results/integration_of_HIC/", "DAR_results_full.xlsx"), overwrite = T)
#openxlsx::saveWorkbook(XLSX_OUT, file.path("~/workspace/neuroblastoma/results/integration_of_HIC/", "DAR_results_only_identified_enh-genes.xlsx"), overwrite = T)



metadata_heatmap <- as.data.frame(colData(ATAC_dds))
heatmap_counts <- SummarizedExperiment::assay(vsd)

rownames(heatmap_counts) <- paste0( vsd@rowRanges$peak_id,"_", vsd@rowRanges$ENH)

# removed xperiment batch effects! ? use experimen+cell_line removed???

#here we need to select peaks that are significant
heatmap_counts <- heatmap_counts[sapply(strsplit(rownames(heatmap_counts), "_"),"[[",1) %in% ATAC_dds_results$results_signif$peak_id,]

annotation_col <- metadata_heatmap %>%
  dplyr::select(cell_line, phenotype) %>% 
  dplyr::arrange(phenotype, cell_line)

heatmap_counts<- heatmap_counts[, match(rownames(annotation_col), colnames(heatmap_counts))]

ensembl2symbol_annot <- ATAC_dds_results$results_all %>%
  dplyr::select(peak_id, gencode_gene_name)

#color.scheme <- rev(RColorBrewer::brewer.pal(8,"RdBu")) # generate the color scheme to use

#Alternative color schemes
#color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012", "#9B2226")
color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
#color.scheme <- c("#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")

ann_colors = list(
  cell_line = c( CLBM = "#005f73", NB10 = "#0a9396", NB8 = "#DD3344", SH = "#FF9F1C"),
  cell_line_id = c(A = "#E9D8A6", M = "#D9BE6D")
)

heatmap <- pheatmap::pheatmap(heatmap_counts,
                              main = "Heatmap of signif. DAR ADRN vs MES.",
                              scale = "row",
                              annotation_col = annotation_col,
                              annotation_colors = ann_colors,
                              #annotation_row = row_annot_symbols,
                              show_colnames = FALSE,
                              show_rownames = T,
                              cluster_cols = FALSE,
                              #cluster_rows = counts_deg_ord_row_cor_hclust,
                              color = color.scheme,
                              #cutree_rows = 5,
                              #gaps_col = gap_cols,
                              #labels_row = make_bold_names(vsd_tp2_vs_tp1_inclEXP9_7_counts_deg_ord_symbols, rownames, rc_names = names(genes_of_interest)),
                              fontsize = 10, fontsize_row = 4) #height=10, cellwidth = 11, cellheight = 11






##### see what genes are enriched with associated peaks:
ATAC_UP_enriched <- ATAC_dds_results$results_signif %>% filter(log2FoldChange > 0 ) %>% select(ENH)
ATAC_UP_enriched <- unlist(strsplit(unlist(ATAC_UP_enriched), ";"))
names(ATAC_UP_enriched) <- NULL

ATAC_DOWN_enriched <- ATAC_dds_results$results_signif %>% filter(log2FoldChange < 0 ) %>% select(ENH)
ATAC_DOWN_enriched <- unlist(strsplit(unlist(ATAC_DOWN_enriched), ";"))
names(ATAC_DOWN_enriched) <- NULL



sort(table(ATAC_UP_enriched), decreasing = T)
sort(table(ATAC_DOWN_enriched), decreasing = T)

intersect(ATAC_DOWN_enriched, ATAC_UP_enriched)

ATAC_UP_enriched_unique <- setdiff(ATAC_UP_enriched, intersect(ATAC_DOWN_enriched, ATAC_UP_enriched))
ATAC_DOWN_enriched_unique <- setdiff(ATAC_DOWN_enriched, intersect(ATAC_DOWN_enriched, ATAC_UP_enriched))


#optionally - check if the genes are expressed in the RNA-seq. For this run RNA-seq script and get the results deg_results

ATAC_UP_enriched_unique <- ATAC_UP_enriched_unique[ATAC_UP_enriched_unique %in% deg_results$results_all$gene_symbol]
ATAC_DOWN_enriched_unique <-ATAC_DOWN_enriched_unique[ATAC_DOWN_enriched_unique %in% deg_results$results_all$gene_symbol]


unique(ATAC_DOWN_enriched)
unique(ATAC_UP_enriched)



gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean=TRUE) 
#gs_C1 <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C1"), clean=TRUE) 
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
#gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean=TRUE) 
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 


C5_GOBP_UP <- hypeR::hypeR(signature = ATAC_UP_enriched_unique, 
                           genesets = gs_C5_GOBP, 
                           test="hypergeometric")
C5_GOBP_plot_UP <- hypeR::hyp_dots(C5_GOBP_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: MES vs ADR YAP/TAZ/TEAD associated genes UP regulated") +theme_bw()

C5_GOBP_DOWN <- hypeR::hypeR(signature = ATAC_DOWN_enriched_unique, 
                             genesets = gs_C5_GOBP, 
                             test="hypergeometric")
C5_GOBP_plot_DOWN <- hypeR::hyp_dots(C5_GOBP_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: MES vs ADR YAP/TAZ/TEAD associated genes DOWN regulated") +theme_bw()


C5_GOCC_UP <- hypeR::hypeR(signature = ATAC_UP_enriched_unique, 
                           genesets = gs_C5_GOCC, 
                           test="hypergeometric")
C5_GOCC_plot_UP <- hypeR::hyp_dots(C5_GOCC_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOCC: UP MES vs ADR YAP/TAZ/TEAD associated genes  ") +theme_bw()


C5_GOCC_DOWN <- hypeR::hypeR(signature = ATAC_DOWN_enriched_unique, 
                             genesets = gs_C5_GOCC, 
                             test="hypergeometric")
C5_GOCC_plot_DOWN <- hypeR::hyp_dots(C5_GOCC_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOCC: DOWN MES vs ADR YAP/TAZ/TEAD associated genes ") +theme_bw()



C5_GOMF_UP <- hypeR::hypeR(signature = ATAC_UP_enriched_unique, 
                           genesets = gs_C5_GOMF, 
                           test="hypergeometric")                        
C5_GOMF_plot_UP <- hypeR::hyp_dots(C5_GOMF_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOMF: UP MES vs ADR YAP/TAZ/TEAD associated genes ") +theme_bw()


C5_GOMF_DOWN <- hypeR::hypeR(signature = ATAC_DOWN_enriched_unique, 
                             genesets = gs_C5_GOMF, 
                             test="hypergeometric")
C5_GOMF_plot_DOWN <- hypeR::hyp_dots(C5_GOMF_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOMF: DOWN MES vs ADR YAP/TAZ/TEAD associated genes ") +theme_bw()



C2_kegg_UP <- hypeR::hypeR(signature = ATAC_UP_enriched_unique, 
                           genesets = gs_C2_kegg, 
                           test="hypergeometric")
C2_kegg_plot_UP <- hypeR::hyp_dots(C2_kegg_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="KEGG:UP MES vs ADR YAP/TAZ/TEAD associated genes") +theme_bw()


C2_kegg_DOWN <- hypeR::hypeR(signature = ATAC_DOWN_enriched_unique, 
                             genesets = gs_C2_kegg, 
                             test="hypergeometric")
C2_kegg_plot_DOWN <- hypeR::hyp_dots(C2_kegg_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="KEGG: DOWN MES vs ADR YAP/TAZ/TEAD associated genes") +theme_bw()




#### Identify TEAD or AP-1 binding sites interactors #####
# 1) Identify TFBSs, 
# 2) subset only the ATAC-seq peaks that contain predicted binding sites
# 3) find potential interactors using Hi-C data.
# 4) subset only the ones that are located in a specific range from TSSs
# 5) check RNA expressoin of these genes


# Identifying TFBSs
# First, we load the PWM sets from JASPAR2022 DB.
# 9606 is a code for human
opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["species"]] <- "9606"
opts[["collection"]] <- "CORE"
opts[["all_versions"]] <- FALSE
motifsToScan <- TFBSTools::getMatrixSet(JASPAR2022, opts)
# names of TEAD motifs are taken from here:
# https://jaspar.genereg.net/search?page=1&version=all&tax_group=all&class=all&family=all&tax_id=all&type=all&collection=all&q=TEAD
#  MA0808.1 - TEAD3
#  MA1121.1 - TEAD2
#  MA0090.3 - TEAD1
#  MA0809.2 - TEAD4
#  MA0099.3 - FOS:JUN
#  MA1128.1 - FOSL1:JUN
motifsToScan <- motifsToScan[names(motifsToScan) %in% c( "MA0090.3", "MA0808.1",  "MA0809.2", "MA1121.1", "MA0099.3", "MA1128.1")]


ATAC_dds <- readRDS(file = "~/workspace/neuroblastoma/RDSs/ATAC_dds.RDS")

ATAC_dds <- estimateSizeFactors(ATAC_dds)
ATAC_dds <- DESeq(ATAC_dds)

vsd <- vst(ATAC_dds, blind = F)

# Now - extracting results

resultsNames(ATAC_dds)
coeff_name <- 'phenotype_M_vs_A'
cond_numerator <-  "A"
cond_denominator <-  "M"
cond_variable <-  "phenotype"
padj_cutoff = 0.05
log2FC_cutoff = 0.58

ATAC_dds_results <- extract_results_DDS(dds_object = ATAC_dds,
                                        coeff_name = coeff_name,
                                        cond_numerator = cond_numerator,
                                        cond_denominator = cond_denominator,
                                        cond_variable = cond_variable,
                                        padj_cutoff = 0.05,
                                        log2FC_cutoff = 0.58)

#keep only DARs
ATAC_dds_signif <- ATAC_dds[ATAC_dds@rowRanges@elementMetadata@listData[["peak_id"]] %in% ATAC_dds_results$results_signif$peak_id,]

# here we check if there are any peaks that have less than 5 reads across all samples.
dim(ATAC_dds_signif)
counts_consensus_filt <- ATAC_dds_signif[rowSums(assay(ATAC_dds_signif)) > 5, ]
dim(counts_consensus_filt)

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

#there are two ways to extract the binding sites. I will use the position information. Save this information as a table and use Julia to find overlaps with the HiCdata.
# motif_matches <- motifmatchr::matchMotifs(pwms = motifsToScan, 
#                                           subject = counts_consensus_filt, 
#                                           genome = BSgenome.Hsapiens.UCSC.hg38, 
#                                           out = "matches")

motif_matches <- motifmatchr::matchMotifs(pwms = motifsToScan, 
                                          subject = counts_consensus_filt, 
                                          genome = BSgenome.Hsapiens.UCSC.hg38, 
                                          out = "positions")


write.csv2(motif_matches$MA0808.1, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites.csv", quote = F, row.names = F) 
write.csv2(motif_matches$MA1121.1, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites.csv", quote = F, row.names = F) 
write.csv2(motif_matches$MA0090.3, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites.csv", quote = F, row.names = F) 
write.csv2(motif_matches$MA0809.2, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites.csv", quote = F, row.names = F) 
write.csv2(motif_matches$MA0099.3, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites.csv", quote = F, row.names = F) 
write.csv2(motif_matches$MA1128.1, file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA1128.1_binding_sites.csv", quote = F, row.names = F) 



# this file was processed using Python script ...../neuroblastoma/03_overlap_between_tead_and_hic_MSC.py
# This gave not amazing results, so we will try other HiC dataset.
# or alternatively - ...../neuroblastoma/03_overlap_between_tead_and_hic_MSC.py
# Now, process output ..../neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results.tsv


#####
#here - to alterate the binding sites for MSC
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites_results.csv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA1128.1_binding_sites_results.csv", header = T)

#here - to alterate the binding sites for IMR90_fibroblasts
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites_results.tsv", header = T)
#TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites_results.csv", header = T)
TF_bindingSites <- read.table("~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA1128.1_binding_sites_results.csv", header = T)



# HIC_bin2 - is the column where target binds. We need to take these intervals and annotate. Bins are 5kb long
#convert it to GRanges object 
# Chr- chromosome name
# HIC_bin2 - this is the location of an interactor of our target. - extracted from HIC
# -no strand
# TArget_start - it's the location of the motif in our ATAC-seq data
# TArget_end - it's the location of the motif in our ATAC-seq data
# HIC_reference_Start - this is the location of the bin in HiC data that corresponds to our target
# HIC_reference_end - this is the location of the bin in HiC data that corresponds to our target (all bins are 5kb)
# HIC_distfoldchange - score from HIC


TF_bindingSites_GR <- GRanges(seqnames = TF_bindingSites$Chr, 
                              ranges = IRanges(TF_bindingSites$HIC_bin2, width = 5000),
                              strand = NULL,
                              Target_Start = TF_bindingSites$Target_Start,
                              Target_End = TF_bindingSites$Target_End,
                              HIC_Reference_Start =  TF_bindingSites$HIC_Reference_Start,
                              HIC_Reference_End =  TF_bindingSites$HIC_Reference_End, 
                              HIC_distfoldchange =  TF_bindingSites$HIC_distfoldchange)

# load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

import::from(.from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"), 
             "filterDatasets",
             "generateResults",
             .character_only=TRUE) 

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
## define promoter regions as +/- 1kb (3kb was used before)
promoter = ChIPseeker::getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
#annotate peaks. promoter determined as ±1kb around TSS
peakAnno  = annotatePeak(TF_bindingSites_GR,tssRegion=c(-1000,1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peak.anno = as.data.frame(peakAnno)
head(peak.anno)

#leave only promoters
peak.anno <- peak.anno[str_detect(peak.anno$annotation, pattern = "Promoter"),]
peak.anno <- peak.anno[!is.na(peak.anno$SYMBOL),]

#write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA0808.1_binding_sites_results_annotated_1kb_prom.tsv")
#write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA1121.1_binding_sites_results_annotated.tsv")
#write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA0090.3_binding_sites_results_annotated.tsv")
#write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA0809.2_binding_sites_results_annotated.tsv")
#write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA0099.3_binding_sites_results_annotated_1kb_prom.tsv")
write.csv2(peak.anno, file = "~/workspace/neuroblastoma/results/integration_of_HIC/IMR90_fibr/MA1128.1_binding_sites_results_annotated_1kb_prom.tsv")
#now - use RNA-seq data to see what genes are expressed.

#load libraries
import::from(.from = DESeq2, .all=TRUE)
import::from(.from = here::here("~/workspace/neuroblastoma/resources/UtilityScriptsRNA-seq.R"), 
             "generateResults",
             "generateEnsemblAnnotation",
             .character_only=TRUE) # used for filtering

#load the RNA-seq data
param_list <- list(
  abs_filt_samples=2,
  padj_cutoff = 0.05,
  log2FC_cutoff = 0.58,
  var_expl_needed = 0.6,
  biomart_host="http://www.ensembl.org", 
  biomart_dataset="hsapiens_gene_ensembl", 
  biomart_Ens_version="Ensembl Genes 109"
)


#loading the DESeq2 file
rnaseq_dds <- readRDS("~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_deseq_data_set.rds")
design(rnaseq_dds) <- as.formula("~ cell_line + cell_type")

#load annotation that was used during pre-processing of the data
annotationData <- read.table(file = "~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_annotation_gene.tsv",
                             sep = "\t",
                             header = TRUE)
#this is necessary to rename columns so that extraction of the data works correctly
colnames(annotationData)[c(5,7)] <- c("ensembl_id", "gene_symbol")

#removeing NB6 sample - this is intermixed sample, the quality is not great
rnaseq_dds <- rnaseq_dds[,rnaseq_dds$cell_line != "STA_NB_6"]
rnaseq_dds$cell_line  <- droplevels(rnaseq_dds$cell_line)

#fillterring lowly expressed genes
rnaseq_dds_filt <- filterDatasets(rnaseq_dds, 
                                  abs_filt = TRUE, 
                                  abs_filt_samples = param_list$abs_filt_samples)

rnaseq_dds_filt <- DESeq2::estimateSizeFactors(rnaseq_dds_filt)
rnaseq_dds_filt <- DESeq2::DESeq(rnaseq_dds_filt)

vsd <- DESeq2::vst(rnaseq_dds_filt, blind = TRUE) # blind = TRUE for QC

deg_results <- generateResults(dds_object = rnaseq_dds_filt, 
                               coeff_name = "cell_type_MES_vs_ADR",
                               cond_numerator = "MES", 
                               cond_denominator = "ADR",
                               cond_variable="cell_type",
                               ensemblAnnot = annotationData)

#filter out genes that are not in the annotated list 
deg_results$results_signif <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(peak.anno$ENSEMBL),]


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
  cell_line = c( CLB_Ma = "#005f73", STA_NB_10 = "#0a9396", STA_NB_8 = "#DD3344", SK_N_SH = "#FF9F1C"),
  cell_line_id = c(ADR = "#E9D8A6", MES = "#D9BE6D")
)

heatmap<- pheatmap::pheatmap(heatmap_counts_deg_ord,
                             main = "Heatmap of signif. DEG",
                             scale = "row",
                             annotation_col = annotation_col,
                             annotation_colors = ann_colors,
                             #annotation_row = row_annot_symbols,
                             show_colnames = FALSE,
                             show_rownames = FALSE,
                             cluster_cols = FALSE,
                             #cluster_rows = counts_deg_ord_row_cor_hclust,
                             color = color.scheme,
                             #cutree_rows = 5,
                             #gaps_col = gap_cols,
                             fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11


########
# now - compare all peaks genes that are affected by TEAD1,2,3,4
########

#TEAD3 <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results_annotated.tsv")
TEAD3 <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0808.1_binding_sites_results_annotated_1kb_prom.tsv")
TEAD2 <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA1121.1_binding_sites_results_annotated.tsv")
TEAD1 <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0090.3_binding_sites_results_annotated.tsv")
TEAD4 <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0809.2_binding_sites_results_annotated.tsv")
FOS_JUN <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA0099.3_binding_sites_results_annotated_1kb_prom.tsv")
FOSL_JUN <- read.csv2( file = "~/workspace/neuroblastoma/results/integration_of_HIC/MA1128.1_binding_sites_results_annotated_1kb_prom.tsv")

rnaseq_dds <- readRDS("~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_deseq_data_set.rds")
design(rnaseq_dds) <- as.formula("~ cell_line + cell_type")

#load annotation that was used during pre-processing of the data
annotationData <- read.table(file = "~/workspace/neuroblastoma/data_soren/ATAC_BSA_0729_KB_NB_plasticity/hg38/rnaseq_deseq_global/rnaseq_deseq_global_annotation_gene.tsv",
                             sep = "\t",
                             header = TRUE)
#this is necessary to rename columns so that extraction of the data works correctly
colnames(annotationData)[c(5,7)] <- c("ensembl_id", "gene_symbol")

#removeing NB6 sample - this is intermixed sample, the quality is not great
rnaseq_dds <- rnaseq_dds[,rnaseq_dds$cell_line != "STA_NB_6"]
rnaseq_dds$cell_line  <- droplevels(rnaseq_dds$cell_line)

#fillterring lowly expressed genes
rnaseq_dds_filt <- filterDatasets(rnaseq_dds, 
                                  abs_filt = TRUE, 
                                  abs_filt_samples = param_list$abs_filt_samples)

rnaseq_dds_filt <- DESeq2::estimateSizeFactors(rnaseq_dds_filt)
rnaseq_dds_filt <- DESeq2::DESeq(rnaseq_dds_filt)

vsd <- DESeq2::vst(rnaseq_dds_filt, blind = TRUE) # blind = TRUE for QC

deg_results <- generateResults(dds_object = rnaseq_dds_filt, 
                               coeff_name = "cell_type_MES_vs_ADR",
                               cond_numerator = "MES", 
                               cond_denominator = "ADR",
                               cond_variable="cell_type",
                               ensemblAnnot = annotationData)

#filter out genes that are not in the annotated list  

deg_TEAD1 <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(TEAD1$ENSEMBL),]
deg_TEAD2 <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(TEAD2$ENSEMBL),]
deg_TEAD3 <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(TEAD3$ENSEMBL),]
deg_TEAD4 <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(TEAD4$ENSEMBL),]
deg_FOS_JUN <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(FOS_JUN$ENSEMBL),]
deg_FOSL_JUN <- deg_results$results_signif[deg_results$results_signif$ensembl_id %in% unique(FOSL_JUN$ENSEMBL),]


TEADS <- list(
  TEAD1=deg_TEAD1$gene_symbol,
  TEAD2=deg_TEAD2$gene_symbol,
  TEAD3=deg_TEAD3$gene_symbol,
  TEAD4=deg_TEAD4$gene_symbol,
  FOSJUN=deg_FOS_JUN$gene_symbol,
  FOSLJUN=deg_FOSL_JUN$gene_symbol
)

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)

ggvenn(
  TEADS, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)


TEADS_intersected <- Reduce(intersect, TEADS)


# just for a small test
#TEADS_intersected <- TEAD3$SYMBOL



deg_results$results_signif <- deg_results$results_signif[deg_results$results_signif$gene_symbol %in% unique(TEADS_intersected),]

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
  cell_line = c( CLB_Ma = "#005f73", STA_NB_10 = "#0a9396", STA_NB_8 = "#DD3344", SK_N_SH = "#FF9F1C"),
  cell_line_id = c(ADR = "#E9D8A6", MES = "#D9BE6D")
)

heatmap<- pheatmap::pheatmap(heatmap_counts_deg_ord,
                             main = "Heatmap of signif. DEG",
                             scale = "row",
                             annotation_col = annotation_col,
                             annotation_colors = ann_colors,
                             #annotation_row = row_annot_symbols,
                             show_colnames = FALSE,
                             show_rownames = FALSE,
                             cluster_cols = FALSE,
                             #cluster_rows = counts_deg_ord_row_cor_hclust,
                             color = color.scheme,
                             #cutree_rows = 5,
                             #gaps_col = gap_cols,
                             #labels_row = make_bold_names(vsd_tp2_vs_tp1_inclEXP9_7_counts_deg_ord_symbols, rownames, rc_names = names(genes_of_interest)),
                             fontsize = 10, fontsize_row = 10) #height=10, cellwidth = 11, cellheight = 11



deg_results <- generateResults(dds_object = rnaseq_dds_filt, 
                               coeff_name = "cell_type_MES_vs_ADR",
                               cond_numerator = "MES", 
                               cond_denominator = "ADR",
                               cond_variable="cell_type",
                               ensemblAnnot = annotationData)


RNA_TARGETS_UP_enriched_unique <- deg_results$results_signif$gene_symbol[ deg_results$results_signif$gene_symbol %in% TEADS_intersected & deg_results$results_signif$log2FoldChange > 0]
RNA_TARGETS_DOWN_enriched_unique <- deg_results$results_signif$gene_symbol[ deg_results$results_signif$gene_symbol %in% TEADS_intersected & deg_results$results_signif$log2FoldChange < 0]

gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean=TRUE) 
#gs_C1 <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C1"), clean=TRUE) 
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean=TRUE) 
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean=TRUE) 
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean=TRUE) 
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean=TRUE) 
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean=TRUE) 


C5_GOBP_UP <- hypeR::hypeR(signature = RNA_TARGETS_UP_enriched_unique, 
                           genesets = gs_C5_GOBP, 
                           test="hypergeometric")
C5_GOBP_plot_UP <- hypeR::hyp_dots(C5_GOBP_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: TEAD associated genes UP regulated") +theme_bw()

C5_GOBP_DOWN <- hypeR::hypeR(signature = RNA_TARGETS_DOWN_enriched_unique, 
                             genesets = gs_C5_GOBP, 
                             test="hypergeometric")
C5_GOBP_plot_DOWN <- hypeR::hyp_dots(C5_GOBP_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOBP: TEAD associated genes DOWN regulated") +theme_bw()


C5_GOCC_UP <- hypeR::hypeR(signature = RNA_TARGETS_UP_enriched_unique, 
                           genesets = gs_C5_GOCC, 
                           test="hypergeometric")
C5_GOCC_plot_UP <- hypeR::hyp_dots(C5_GOCC_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOCC: UP TEAD associated genes  ") +theme_bw()


C5_GOCC_DOWN <- hypeR::hypeR(signature = RNA_TARGETS_DOWN_enriched_unique, 
                             genesets = gs_C5_GOCC, 
                             test="hypergeometric")
C5_GOCC_plot_DOWN <- hypeR::hyp_dots(C5_GOCC_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOCC: DOWN TEAD associated genes ") +theme_bw()


C5_GOMF_UP <- hypeR::hypeR(signature = RNA_TARGETS_UP_enriched_unique, 
                           genesets = gs_C5_GOMF, 
                           test="hypergeometric")                        
C5_GOMF_plot_UP <- hypeR::hyp_dots(C5_GOMF_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOMF: UP TEAD associated genes ") +theme_bw()


C5_GOMF_DOWN <- hypeR::hypeR(signature = RNA_TARGETS_DOWN_enriched_unique, 
                             genesets = gs_C5_GOMF, 
                             test="hypergeometric")
C5_GOMF_plot_DOWN <- hypeR::hyp_dots(C5_GOMF_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="GOMF: DOWN TEAD associated genes ") +theme_bw()



C2_kegg_UP <- hypeR::hypeR(signature = RNA_TARGETS_UP_enriched_unique, 
                           genesets = gs_C2_kegg, 
                           test="hypergeometric")
C2_kegg_plot_UP <- hypeR::hyp_dots(C2_kegg_UP, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="KEGG:UP TEAD associated genes") +theme_bw()


C2_kegg_DOWN <- hypeR::hypeR(signature = RNA_TARGETS_DOWN_enriched_unique, 
                             genesets = gs_C2_kegg, 
                             test="hypergeometric")
C2_kegg_plot_DOWN <- hypeR::hyp_dots(C2_kegg_DOWN, merge=TRUE, fdr=0.05, top = 20, abrv=70, val="fdr", title="KEGG: DOWN TEAD associated genes") +theme_bw()

