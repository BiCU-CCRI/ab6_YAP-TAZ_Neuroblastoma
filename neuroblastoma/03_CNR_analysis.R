###
# Title: C&R-seq analysis for Sören (Kaan's group)
# Author: Aleksandr Bykov
#

# fix for diffBind Greylist filterring. It doesn't work if the original file containing scaffold chromosomes.
pv.countGreylistEdited <- function(bamfile, pv, ktype) {
  # gl <- new("GreyList", karyotype = ktype[pv$chrmap, ]) #previous line
  # edit to restrict to just chromosomes included in the ktype object
  gl <- new("GreyList", karyotype = ktype[intersect(pv$chrmap, names(ktype)), ])
  gl <- GreyListChIP::countReads(gl, bamfile)
  return(gl)
}
environment(pv.countGreylistEdited) <- asNamespace("DiffBind")
assignInNamespace("pv.countGreylist", pv.countGreylistEdited, ns = "DiffBind")

# set up the environment
res_dir <- "/home/rstudio/workspace/neuroblastoma/results/CnR/"

#### Loading libraries ####
import::from(
  .from = "~/workspace/neuroblastoma/resources/utilityScripts.R",
  "generatePCA",
  "extract_results_DDS",
  "meanExprsPerGroup",
  "extract_results_DDS_HIC",
  "gseaplot3",
  "tableGrob2",
  "gsInfo",
  "chipEnrichAndExport"
)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(stringr)
library(chipenrich)
library(DiffBind)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(biomaRt)

library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg38)
library(genio)
library(org.Hs.eg.db)
library(viridis)

#library(JASPAR2022) # fix for JASPAR library
system("mkdir /home/rstudio/.cache/R/BiocFileCache/")
download.file(url = "https://jaspar2022.genereg.net/download/database/JASPAR2022.sqlite",
              destfile = "/home/rstudio/.cache/R/BiocFileCache/JASPAR2022.sqlite")
JASPAR2022 <-  "/home/rstudio/.cache/R/BiocFileCache/JASPAR2022.sqlite"

# TODO fix the timeout
mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)

import::from(.from = magrittr, "%>%")
import::from(.from = IRanges, IRanges)
import::from(.from = GenomicRanges, GRanges)
import::from(.from = rtracklayer, liftOver, import.chain)
import::from(.from = openxlsx, addWorksheet, writeData, saveWorkbook, createWorkbook)
import::from(.from = openxlsx2, read_xlsx)
import::from(.from = TFBSTools, getMatrixSet)
import::from(.from = chromVAR, addGCBias)
import::from(.from = RColorBrewer, brewer.pal)
import:::here(.from = DOSE, gseaScores)


futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# loading the GO data
gs_hallmark    <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"),                               clean = TRUE)
gs_C2_kegg     <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG",     clean = TRUE)
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
gs_C5_GOBP     <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP",       clean = TRUE)
gs_C5_GOCC     <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC",       clean = TRUE)
gs_C5_GOMF     <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF",       clean = TRUE)

system("zip -F neuroblastoma/data/CnR/consensus_peaks.zip --out neuroblastoma/data/CnR/consensus_peaks_u.zip")
system("unzip neuroblastoma/data/CnR/consensus_peaks_u.zip")
system("rm neuroblastoma/data/CnR/consensus_peaks_u.zip")

data_folder <- "/home/rstudio/workspace/neuroblastoma/data/CnR"

files <- unlist(list.files(path = file.path(data_folder, "consensus_peaks"), pattern = "dba", full.names = TRUE))
sample_names <- stringr::str_extract(string = files, 
                                     pattern = ".*(dba_)(.*).RData", group = 2)
DBobj_list <- list()
DBobj_list_cons <- list()
for (target_prot in unique(sample_names)) {
  print(paste("File", target_prot, "exists, loading ..."))
  fl_path <- file.path(data_folder, paste0("consensus_peaks/dba_", target_prot))
  dfbobj <- dba.load(fl_path, dir = "", pre = "")
  DBobj_list[[target_prot]] <- dfbobj
  DBobj_list_cons[[target_prot]] <- dba.peakset(dfbobj, bRetrieve = TRUE)
}



# ######
# dba.overlap(DBobj_list$H3K27ac, mode = DBA_OLAP_RATE)
# dba.overlap(DBobj_list$H3K4me1, mode = DBA_OLAP_RATE)
# dba.overlap(DBobj_list$Jun, mode = DBA_OLAP_RATE)
# #####

# # Load all peaks from .narrowPeaks files
# folder_path <- "~/workspace/neuroblastoma/data/CnR/peaks/"
# files <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.narrowPeak$", full.names = TRUE))
# summits <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.macs2_summits.bed$", full.names = TRUE))
# 
# # assemble the samples object
# samples <- data.frame(list(SampleID = str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R(\\d)")))
# samples$Tissue <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)")
# samples$Factor <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R", group = 3)
# samples$Condition <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)", group = 2)
# samples$Replicate <- str_extract(files, pattern = "_R(\\d)", group = 1)
# samples$Peaks <- files
# samples$PeakCaller <- "narrow"
# samples$PeakSummits <- summits
# peaks_count <- c()
# peaks_count <- sapply(files, function(x) genio::count_lines(x) - 1)
# samples$Peaks_count <- peaks_count
# 
# # Loading all narrowPeaks
# DBobj_list_cons <- list()
# for (target_prot in unique(samples$Factor)) {
#   print(c("Current protein is ", target_prot))
#   dfbobj <- dba(sampleSheet = samples %>% filter(., Factor == target_prot))
#   DBobj_list_cons[[target_prot]] <- dba.peakset(dfbobj, bRetrieve = TRUE)
# }

peakAnnoList <- list()
# check features distribution
for (factor in names(DBobj_list_cons)){
  test_region <- DBobj_list_cons[[factor]]
  seqlevelsStyle(test_region) <- "UCSC"
  peakAnnoList[[factor]] <- ChIPseeker::annotatePeak(peak = test_region, 
                                         tssRegion = c(-3000, 3000),
                                         TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                         annoDb = "org.Hs.eg.db")
}
plotAnnoBar(peakAnnoList)

# technical analysis - Check overlap between YAP TAZ and Jun peaks
TAZ_Peaks <- DBobj_list_cons$TAZ
YAP_Peaks <- DBobj_list_cons$YAP
Jun_Peaks <- DBobj_list_cons$Jun

ol <- ChIPpeakAnno::findOverlapsOfPeaks(TAZ_Peaks, YAP_Peaks)
ChIPpeakAnno::makeVennDiagram(ol,
  fill = c("#f5c9b1", "#f1b5ab"), # circle fill color
  col = c("black", "black"), # circle border color
  cat.col = c("#D55E00", "#0072B2")
)

# correlation of signal
for_corr_plot <- ol$overlappingPeaks$`TAZ_Peaks///YAP_Peaks` %>%
  dplyr::select(
    CLB.Ma.M_TAZ_R1, CLB.Ma.M_TAZ_R2, SK.N.SH.M_TAZ_R1, SK.N.SH.M_TAZ_R2,
    CLB.Ma.M_YAP_R1, CLB.Ma.M_YAP_R2, SK.N.SH.M_YAP_R1
  )

for_corr_plot <- data.frame(
  TAZ = apply(for_corr_plot[, 1:4], 1, mean),
  YAP = apply(for_corr_plot[, 5:7], 1, mean)
)

# Plot correlation for TAZ and YAP signal
ggplot(for_corr_plot) +
  aes(x = TAZ, y = YAP) +
  geom_point(
    shape = "circle", size = 1.5,
    colour = "#112446"
  ) +
  theme_minimal() +
  ggpubr::stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01) +
  geom_smooth(method = "lm")


ol <- ChIPpeakAnno::findOverlapsOfPeaks(ol$peaklist$`TAZ_Peaks///YAP_Peaks`, Jun_Peaks)
ChIPpeakAnno::makeVennDiagram(ol,
  fill = c("#f1b5ab", "#948bbe"), # circle fill color
  col = c("black", "black"), # circle border color
  cat.col = c("#0072B2", "black")
)

ol <- ChIPpeakAnno::findOverlapsOfPeaks(TAZ_Peaks, YAP_Peaks, Jun_Peaks)
ChIPpeakAnno::makeVennDiagram(ol)

##### Running Gene set enrichment using the chipenrich package ##########
# turn RNA-seq results to Terms table
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = mes_adrn_gene_list$ensembl_id,
  mart = mart
)
mes_adrn_gene_list <-
  merge.data.frame(
    x = mes_adrn_gene_list,
    y = genes,
    by.x = "ensembl_id",
    by.y = "ensembl_gene_id"
  ) %>%
  dplyr::select(Term, entrezgene_id) %>%
  dplyr::rename(gs_id = Term, gene_id = entrezgene_id) %>%
  dplyr::filter(!is.na(gene_id))

write.table(
  x = mes_adrn_gene_list,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)

locusdef <-  "5kb"

#our_rna_seq_terms <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv"
#our_rna_seq_terms <- mes_adrn_gene_list
our_rna_seq_terms <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv"

for(TF_name in names(DBobj_list)){
  #TF_name <- "YAP"  
  print(paste0("Processing ", TF_name))
  
  # subset peaks that up-regulated in MES samples
  p_MES_up <- as.data.frame(dba.report(DBobj_list[[TF_name]], contrast = 2)) %>%
    dplyr::filter(Fold > 0) %>%
    dplyr::select(seqnames, start, end)
  
  # process up-regulated peaks
  # remove scaffolds - their name is longer than 2 charcaters
  p_MES_up <- p_MES_up %>% dplyr::filter(nchar(as.character(seqnames)) <= 2)
  p_MES_up$seqnames <- droplevels(p_MES_up$seqnames)
  p_MES_up <- p_MES_up %>% dplyr::mutate(seqnames = paste0("chr", seqnames))
  colnames(p_MES_up)[1] <- "chrom"
  #plot_dist_to_tss(peaks = p_MES_up, genome = "hg38")
  
  # subset peaks that down-regulated in MES samples
  p_MES_dwn <- as.data.frame(dba.report(DBobj_list[[TF_name]], contrast = 2)) %>%
    filter(Fold < 0) %>%
    dplyr::select(seqnames, start, end)
  
  # process down-regulated peaks
  # remove scaffolds - their name is longer than 2 charcaters
  # add chr tag to seqnames
  
  p_MES_dwn <- p_MES_dwn %>% dplyr::filter(nchar(as.character(seqnames)) <= 2)
  p_MES_dwn$seqnames <- droplevels(p_MES_dwn$seqnames)
  p_MES_dwn <- p_MES_dwn %>% dplyr::mutate(seqnames = paste0("chr", seqnames))
  colnames(p_MES_dwn)[1] <- "chrom"
  #plot_dist_to_tss(peaks = p_MES_dwn, genome = "hg38")
  
  if (nrow(p_MES_up) > 0) {
    chipEnrichAndExport(peaks = p_MES_up,
                        peaksName = "MES", 
                        TF_name = TF_name, 
                        res_dir = res_dir, 
                        genesets = our_rna_seq_terms, 
                        genesets_name = "Our_data_",
                        locusdef = locusdef
                        )
  }
  
  if (nrow(p_MES_dwn) > 0){
    
    chipEnrichAndExport(peaks = p_MES_dwn,
                        peaksName = "ADRN", 
                        TF_name = TF_name, 
                        res_dir = res_dir, 
                        genesets = our_rna_seq_terms, 
                        genesets_name = "Our_data_",
                        locusdef = locusdef
    )
  }
}

# alternative enrichment analysis
####### test how to visualize #########
enrich_results_files <- unlist(list.files(path = res_dir, pattern = "^Enricher_.*.xlsx", full.names = TRUE))
summary_table <- data.frame(Protein = NULL,
                            Data_source = NULL,
                            Selected_peaks = NULL,
                            Description = NULL,
                            P.value = NULL,
                            FDR = NULL,
                            Effect = NULL,
                            Status = NULL,
                            Gene_set_size = NULL,
                            Peaks_in_set = NULL,
                            Odds_ratio = NULL)

#pattern <-  ".*(H3K27ac|H3K4me1|Jun|TAZ|YAP)_(Our_data|Groen).*(ADRN|MES).*"
pattern <-  ".*(H3K27ac|H3K4me1|Jun|TAZ|YAP)_(Our_data|Groen).*(ADRN|MES).*"

for (file_name in enrich_results_files){
  tmp <- openxlsx2::wb_to_df(file_name, sheet = 1)
  
  basename_tmp <- basename(file_name)
  basename_tmp <- str_extract(basename_tmp, 
                              pattern, 
                              group = c(1,2,3))
  
  summary_table_tmp <- data.frame(Protein = basename_tmp[1],
                                  Data_source = basename_tmp[2],
                                  Selected_peaks = basename_tmp[3],
                                  Description = tmp$Description,
                                  P.value = tmp$P.value,
                                  FDR = tmp$FDR,
                                  Effect = tmp$Effect,
                                  Status = tmp$Status,
                                  Gene_set_size = tmp$N.Geneset.Genes,
                                  Peaks_in_set = tmp$N.Geneset.Peak.Genes,
                                  Odds_ratio = tmp$Odds.Ratio)
  
  summary_table <- rbind(summary_table,
                         summary_table_tmp)
  
}

# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_data", Selected_peaks == "MES")
custom_colors <- c("red", colorRampPalette(brewer.pal(7, "Greys"))(100))
custom_breaks <- c(seq(0, 1, length.out = 2), seq(1, 60, length.out = 256))
# making lolipop plot with enrichments
summary_table_MES %>%
  ggplot() +
  geom_bar(position = position_dodge(0.5), 
           width = 0.1, 
           aes(y = Protein, 
               fill = Description, 
               weight = Odds_ratio)) +
  scale_fill_manual(values =  c("navy", "firebrick3")) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(y = Protein, 
                 x = Odds_ratio, 
                 color = Description,
                 size = Peaks_in_set/Gene_set_size*100,
                 fill = -log10(FDR)
                 ), 
             shape = "circle filled",
             position = position_dodge2(0.5)
             ) +
  scale_color_manual(values =  c("navy", "firebrick3")) +
  scale_fill_gradientn(colors = custom_colors, 
                       values = scales::rescale(custom_breaks),
                       limits = c(0, 70)) +
#  scale_fill_viridis(option="viridis")+
  theme_minimal() +
  labs(color = "selected peaks", 
       fill = "-log10(FDR)\n red - non significant (FDR < 0.05)",
       size =  "Percentage of genes in set\n overlaping with gene-set collection ",
       y = "Proteins",
       x = "Odds ratio",
       title = "Enrichment of peaks in MES samples in MES/ADR-specific regions determined from RNA-seq")

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_data", Selected_peaks == "ADRN")
custom_colors <- c("red", colorRampPalette(brewer.pal(7, "Greys"))(100))
custom_breaks <- c(seq(0, 1, length.out = 2), seq(1, 60, length.out = 256))
# making lolipop plot with enrichments
summary_table_ADRN %>%
  ggplot() +
  geom_bar(position = position_dodge(0.5), 
           width = 0.1, 
           aes(y = Protein, 
               fill = Description, 
               weight = Odds_ratio)) +
  scale_fill_manual(values =  c("navy", "firebrick3")) +
  scale_x_continuous(limits = c(0, 15)) +
  ggnewscale::new_scale_fill() +
  geom_point(aes(y = Protein, 
                 x = Odds_ratio, 
                 color = Description,
                 size = Peaks_in_set/Gene_set_size*100,
                 fill = -log10(FDR)
  ), 
  shape = "circle filled",
  position = position_dodge2(0.5)
  ) +
  scale_color_manual(values =  c("navy", "firebrick3")) +
  scale_fill_gradientn(colors = custom_colors, 
                       values = scales::rescale(custom_breaks),
                       limits = c(0, 70)) +
  #  scale_fill_viridis(option="viridis")+
  theme_minimal() +
  labs(color = "selected peaks", 
       fill = "-log10(FDR)\n red - non significant (FDR < 0.05)",
       size =  "Percentage of genes in set\n overlaping with gene-set collection ",
       y = "Proteins",
       x = "Odds ratio",
       title = "Enrichment of peaks in ADRN samples in MES/ADR-specific regions determined from RNA-seq")


# making lolipop plot with enrichments

# 
# summary_table <- summary_table %>% filter(Description == "MES" & FDR < 0.05) %>% 
#   mutate(Protein_Data_source = paste0(Protein, "_", Data_source))
# 
# library(viridis)


# #plot our data
# ggplot(filter(summary_table, Data_source == "Our_data")) +
#   geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf),
#             fill = "blue", alpha = 0.01, inherit.aes = F ) +
#   geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf),
#             fill = "red", alpha = 0.01, inherit.aes = F ) +
#   aes(
#     fill = FDR,
#     size = Peaks_in_set
#   ) +
#   geom_point(aes( x = Effect, y = Protein, fill = FDR), shape = 21, colour = "black", ) +
#   scale_fill_viridis(option = "inferno", direction = 1) +
#   scale_x_continuous(labels = function(x) abs(as.numeric(x)),
#                      name = "      Enrichment in ADRN <--   --> Enrichment in MES") +
#   theme_minimal()

# # plot GROENINGEN data
# ggplot(filter(summary_table, Data_source == "Groen")) +
#   geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf), 
#             fill = "blue", alpha = 0.01, inherit.aes = F ) +
#   geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf), 
#             fill = "red", alpha = 0.01, inherit.aes = F ) +
#   aes(
#     fill = FDR,
#     size = Peaks_in_set
#   ) +
#   geom_point(aes( x = Effect, y = Protein, fill = FDR), shape = 21, colour = "black", ) +
#   scale_fill_viridis(option = "inferno", direction = 1) +
#   scale_x_continuous(labels = function(x) abs(as.numeric(x)),
#                      name = "      Enrichment in ADRN <--   --> Enrichment in MES") +
#   theme_minimal() +
#   ggtitle("GROENINGEN DATA")

# # Split plots
# pdf(file = file.path(res_dir, "/TFs_Hmods_MES_ADRN_our_data.pdf"))
# 
# summary_table %>%
#   filter(Status %in% "enriched") %>%
#   filter(Protein_Data_source %in% c("H3K4me1_Our_data", 
#                                     "H3K27ac_Our_data")) %>%
#   ggplot() +
#   aes(x = Protein, y = Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Histone Modification", title = "Histone Modifications, MES, our Data") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# summary_table %>%
#   filter(Status %in% "enriched") %>%
#   filter(Protein_Data_source %in% c("Jun_Our_data", 
#                                     "TAZ_Our_data",
#                                     "YAP_Our_data")) %>%
#   ggplot() +
#   aes(x = Protein, y = Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Transcription Factors", title = "Transcription Factors, MES, our Data") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# summary_table %>%
#   filter(Status %in% "depleted") %>%
#   filter(Protein_Data_source %in% c("H3K4me1_Our_data", 
#                                     "H3K27ac_Our_data")) %>%
#   ggplot() +
#   aes(x = Protein, y = -Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Histone Modification", title = "Histone Modifications, ADRN, our Data") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# 
# summary_table %>%
#   filter(Status %in% "depleted") %>%
#   filter(Protein_Data_source %in% c("Jun_Our_data", 
#                                     "TAZ_Our_data",
#                                     "YAP_Our_data")) %>%
#   ggplot() +
#   aes(x = Protein, y = -Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Transcription Factors", title = "Transcription Factors, ADRN, our Data") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal()
# 
# dev.off()

# # THE SAME, but for groeningen data
# pdf(file = file.path(res_dir, "plots/TFs_Hmods_MES_ADRN_groeningen.pdf"))
# 
# summary_table %>%
#   filter(Status %in% "enriched") %>%
#   filter(Protein_Data_source %in% c("H3K4me1_Groen", 
#                                     "H3K27ac_Groen")) %>%
#   ggplot() +
#   aes(x = Protein, y = Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Histone Modification", title = "Histone Modifications, MES, Groen") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# summary_table %>%
#   filter(Status %in% "enriched") %>%
#   filter(Protein_Data_source %in% c("Jun_Groen", 
#                                     "TAZ_Groen",
#                                     "YAP_Groen")) %>%
#   ggplot() +
#   aes(x = Protein, y = Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Transcription Factors", title = "Transcription Factors, MES, Groen") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# summary_table %>%
#   filter(Status %in% "depleted") %>%
#   filter(Protein_Data_source %in% c("H3K4me1_Groen", 
#                                     "H3K27ac_Groen")) %>%
#   ggplot() +
#   aes(x = Protein, y = -Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Histone Modification", title = "Histone Modifications, ADRN, Groen") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal() 
# 
# 
# summary_table %>%
#   filter(Status %in% "depleted") %>%
#   filter(Protein_Data_source %in% c("Jun_Groen", 
#                                     "TAZ_Groen",
#                                     "YAP_Groen")) %>%
#   ggplot() +
#   aes(x = Protein, y = -Effect, colour = FDR, size = Peaks_in_set) +
#   ylim(0, 2.5) +
#   labs(y= "Effect", x = "Transcription Factors", title = "Transcription Factors, ADRN, Groen") +
#   geom_point(shape = "circle") +
#   scale_color_gradient() +
#   theme_minimal()
# 
# dev.off()


##############################################
### Check overlapped genes from 3 datasets - YAP CnR, ATAC-seq, RNA-seq (upregul in MES data)
### See how do they compare to Groningen Data
# for YAP, TAZ and Jun
# TODO - label peaks all H3K27ac as ACTIVE


# # TOY TEST:
# peaks1 <- GRanges(seqnames=c(6,6,6,6,6,6),
#                   IRanges(start=c(1,3,5,8,11,14),
#                           end=c(2,4,7,10,13,15),
#                           names=c("p1","p2","p3","p4","p5", "p6")),
#                   strand="+")
# peaks2 <- GRanges(seqnames=c(6,6,6,6,6,6),
#                   IRanges(start=c(6,8,11,14,17,19),
#                           end=c(9,10,12,16,18,20),
#                           names=c("f1","f2","f3","f4","f5", "f6")),
#                   strand="+")
# ChIPpeakAnno::findOverlapsOfPeaks(peaks1, peaks2)

# TODO - label peaks all H3K27ac as ACTIVE

TAZ_Peaks <- dba.peakset(DBobj_list$TAZ, consensus = TRUE, bRetrieve = TRUE)
YAP_Peaks <- dba.peakset(DBobj_list$YAP, consensus = TRUE, bRetrieve = TRUE)
Jun_Peaks <- dba.peakset(DBobj_list$Jun, consensus = TRUE, bRetrieve = TRUE)
H3k27ac_Peaks <- dba.peakset(DBobj_list$H3K27ac, consensus = TRUE, bRetrieve = TRUE)
H3k4me1_Peaks <- dba.peakset(DBobj_list$H3K4me1, consensus = TRUE, bRetrieve = TRUE)

TAZ_Peaks <- DBobj_list_cons$TAZ
YAP_Peaks <- DBobj_list_cons$YAP
Jun_Peaks <- DBobj_list_cons$Jun
H3k27ac_Peaks <- DBobj_list_cons$H3K27ac
H3k4me1_Peaks <- DBobj_list_cons$H3K4me1

ol <- ChIPpeakAnno::findOverlapsOfPeaks(TAZ_Peaks, YAP_Peaks)
YAP_TAZ_peaks <- ol$mergedPeaks
ol <- ChIPpeakAnno::findOverlapsOfPeaks(YAP_TAZ_peaks, Jun_Peaks)
YAP_TAZ_JUN_peaks <- ol$mergedPeaks

YAP_TAZ_JUN_peaks_list <- list(
  Jun_Peaks = Jun_Peaks,
  YAP_TAZ_peaks = YAP_TAZ_peaks,
  YAP_TAZ_JUN_peaks = YAP_TAZ_JUN_peaks
)



promoter <- ChIPseeker::getPromoters(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 3000, downstream = 3000)
promoter <- as.data.frame(promoter)[, 1:3]
promoter$seqnames <- str_replace_all(promoter$seqnames, pattern = "chr", replacement = "")
promoter <- promoter %>% dplyr::arrange(seqnames)
promoter <- promoter[sapply(promoter$seqnames, nchar) <= 2, ]
promoter <- promoter %>%
  mutate(merged = paste0(seqnames, ":", start, ":", end)) %>%
  distinct(merged, .keep_all = T)
promoter <- GRanges(promoter)

tmp <- ChIPpeakAnno::findOverlapsOfPeaks(promoter, H3k27ac_Peaks, H3k4me1_Peaks)
ChIPpeakAnno::makeVennDiagram(tmp)

active_promoters <- tmp$peaklist$`promoter///H3k27ac_Peaks`
poised_promoters <- tmp$peaklist$promoter
active_enhancer  <- tmp$peaklist$`H3k27ac_Peaks///H3k4me1_Peaks`
poised_enhancer  <- tmp$peaklist$H3k4me1_Peaks

piechart <- tibble()
for(peak_names in names(YAP_TAZ_JUN_peaks_list)){
  #peak_names <- "Jun_Peaks"
  peaks <- YAP_TAZ_JUN_peaks_list[[peak_names]]
  
  ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks, active_promoters)
  active_promoters_percent <- length(peaks) - length(ol$peaklist$peaks)
  
  ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks, poised_promoters)
  poised_promoters_percent <- length(peaks) - length(ol$peaklist$peaks)
  
  ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks,  active_enhancer)
  active_enhancer_percent <- length(peaks) - length(ol$peaklist$peaks)
  
  ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks, poised_enhancer)
  poised_enhancer_percent <- length(peaks) - length(ol$peaklist$peaks)
  
  unchar <- length(peaks) - active_promoters_percent - poised_promoters_percent - active_enhancer_percent - poised_enhancer_percent
  
  piechart <- rbind(piechart,
                    tibble(
                      group = c(peak_names),
                      Peaks = c("Unclassified", 
                                "Active Promoters H3K27Ac(+) H3K4me1(-)", 
                                "Inactive/Bivalient Promoters",
                                "Active Enhancers H3K27Ac(+) H3K4me1(+)", 
                                "Poised enhancers H3K27Ac(-) H3K4me1(+)"),
                      number_of_peaks = c(
                        unchar,
                        active_promoters_percent,
                        poised_promoters_percent,
                        active_enhancer_percent,
                        poised_enhancer_percent
                      )
                    )
  )
}

piechart <- piechart %>% group_by(group) %>% mutate(percent = round(number_of_peaks / sum(number_of_peaks) * 100))
  #piechart$percent <- round(piechart$number_of_peaks / sum(piechart$number_of_peaks) * 100, 2)
  
  ggplot(piechart) +
    aes(x = group, y = percent, fill = Peaks) +
    geom_col() +
    scale_fill_brewer(palette = "Accent", direction = 1) +
    theme_minimal() +
    geom_text(aes(label = paste0(percent, "%")),
              position = position_stack(vjust = 0.5), size = 6
    )
  


# Checking the distribution of peaks in gene parts (promoters, 3' UTRs, etc.)
# for (TF in names(DBobj_list)) {
#   tmp <- DBobj_list[[TF]][["DESeq2"]][["DEdata"]]@rowRanges@elementMetadata@listData[["Annotation"]]
#   tmp <- str_replace(tmp, pattern = " (\\(.*\\))", replacement = "")
#   tmp <- str_replace(tmp, pattern = " promoter-TSS..", replacement = "promoter-TSS")
#   tmp <- str_replace(tmp, pattern = "\\.\\d", "")
#   tmp <- tmp[!is.na(tmp)]
# 
#   pie_vector <- c()
# 
#   for (i in unique(tmp)) {
#     pie_vector <- c(pie_vector, sum(i == tmp) / length(tmp))
#   }
#   names(pie_vector) <- unique(tmp)
#   pie(pie_vector)
# }



# now, analyze the diff binding for every protein
for (target_prot in names(DBobj_list)) {
  
  dfbobj <- DBobj_list[[target_prot]]

  hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
  dba.plotHeatmap(dfbobj,
    contrast = 2, correlations = FALSE,
    scale = "row", colScheme = hmap, th = 0.01, bUsePval = F
  )
  plot(dfbobj, contrast = 2)

  rep <- dba.report(dfbobj, contrast = 2)
  repUP <- rep[rep$Fold > 1, ]
  repDWN <- rep[rep$Fold < -1, ]
  repUP <- repUP[order(repUP$Fold, decreasing = TRUE), ]
  repDWN <- repDWN[order(repDWN$Fold, decreasing = FALSE), ]

  write.table(
    as.data.frame(repUP)[, 1:3] %>%
      dplyr::select(seqnames, start, end) %>%
      mutate(seqnames = paste0("chr", seqnames), chain = "."),
    file = paste0("neuroblastoma/results/CnR/", target_prot, "_consenus_peaks_UP.bed"),
    sep = "\t",
    row.names = T,
    col.names = F,
    quote = F
  )

  write.table(
    as.data.frame(repDWN)[, 1:3] %>%
      dplyr::select(seqnames, start, end) %>%
      mutate(seqnames = paste0("chr", seqnames), chain = "."),
    file = paste0("neuroblastoma/results/CnR/", target_prot, "_consenus_peaks_DWN.bed"),
    sep = "\t",
    row.names = T,
    col.names = F,
    quote = F
  )

  if (length(repUP) != 0 & (length(repDWN) != 0)) {
    repList <- GRangesList(UP = repUP, DWN = repDWN)
  } else {
    if (length(repUP) == 0) {
      repList <- GRangesList(DWN = repDWN)
    } else {
      if (length(repDWN) == 0) {
        repList <- GRangesList(UP = repUP)
      }
    }
  }
}


#############
# Peak profiles
# Fix the BAM paths in the original DBobj

pdf(file = file.path(res_dir, "/h3k27ac-h3k4me1-jun-taz-yap.pdf"), width = 12, height = 6)
for (dfbobj in DBobj_list) {
  str_to_be_replaced <- "/home/rstudio/workspace/neuroblastoma/data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/"
  str_replacement <- "/home/rstudio/workspace/neuroblastoma/data/CnR/BAMs/"
  dfbobj$class["bamRead",] <- dfbobj$class["bamRead",] %>% 
    str_replace_all(string = .,
                    pattern = str_to_be_replaced, 
                    replacement = str_replacement)
  dfbobj$class["bamControl",] <- dfbobj$class["bamControl",] %>% 
    str_replace_all(string = .,
                    pattern = str_to_be_replaced, 
                    replacement = str_replacement)
  # str_to_be_replaced <- "/home/rstudio/workspace/neuroblastoma/data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/04_called_peaks/macs2/"
  # str_replacement <- "/home/rstudio/workspace/neuroblastoma/data/CnR/peaks/"
  # dfbobj$samples$Peaks <- dfbobj$samples$Peaks %>% 
  #   str_replace_all(string = .,
  #                   pattern = str_to_be_replaced, 
  #                   replacement = str_replacement)
  #dfbobj <- DBobj_list[[1]]
  rep <- dba.report(dfbobj, contrast = 2)
  
  repUP <- rep[rep$Fold > 2, ]
  repDWN <- rep[rep$Fold < -1, ]
  
  repUP <- repUP[order(repUP$Fold, decreasing = TRUE), ]
  repDWN <- repDWN[order(repDWN$Fold, decreasing = FALSE), ]
  
  repList <- GRangesList(
    UP = repUP,
    DWN = repDWN
  )
  
  rep <- rep[abs(rep$Fold) > 2, ]
  
  if (length(repList$DWN) !=0){
    profiles <- dba.plotProfile(dfbobj,
                                merge = c(DBA_TISSUE, DBA_REPLICATE),
                                contrast = 2,
                                sites = repList)
  }else{
    profiles <- dba.plotProfile(dfbobj,
                                merge = c(DBA_TISSUE, DBA_REPLICATE),
                                contrast = 2,
                                sites = rep)
  }
  print(paste("this is ", unique(dfbobj[["samples"]][["Factor"]])))
  dba.plotProfile(profiles)
}
dev.off()






  
##### Playground
  
# Check the overlap between Groningen data and our data. Check how CnR peaks are overlapping with genes that are
# known to be specific for ADRN or MES
# Load data from Groeningen paper https://www.nature.com/articles/ng.3899 - genes
# that are specific to MES and ADRN cells and see how CnR peaks are placed there
ADR_MES_genes <- read.table("~/workspace/neuroblastoma/resources/mes_adrn_genes.tsv", sep = "\t", header = F, col.names = c("gene", "type"))

# load sets for GO
gs_hallmark <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("H"), clean = TRUE)
gs_C2_kegg <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:KEGG", clean = TRUE)
gs_C2_reactome <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C2"), subcategory = "CP:REACTOME", clean = TRUE)
gs_C5_GOBP <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:BP", clean = TRUE)
gs_C5_GOCC <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:CC", clean = TRUE)
gs_C5_GOMF <- hypeR::msigdb_gsets(species = "Homo sapiens", category = c("C5"), subcategory = "GO:MF", clean = TRUE)

for (DBobj in DBobj_list) {
  # experiment name
  exp_name <- DBobj[["samples"]][["Factor"]][1]
  # initialize pdf document
  pdf(paste0("~/workspace/neuroblastoma/results/20230919/", exp_name, "_results.pdf"))
  # Intersect for MES
  # Extract gene names from the YAP dataset
  Gene_Names <- DBobj[["DESeq2"]][["DEdata"]]@rowRanges@elementMetadata@listData[["Gene.Name"]]
  # Extract indexes of gene names that correspond to a specific criterion. FDR is already < 0.05
  Gene_ids_positive <- as.data.frame(dba.report(DBobj, contrast = 2)) %>% filter(Fold > 0)
  Gene_ids_negative <- as.data.frame(dba.report(DBobj, contrast = 2)) %>% filter(Fold < 0)
  # get the list of gene names
  Gene_Names_MES_positive <- Gene_Names[as.numeric(row.names.data.frame(Gene_ids_positive))]
  Gene_Names_MES_negative <- Gene_Names[as.numeric(row.names.data.frame(Gene_ids_negative))]
  # What are the genes that overlap between Groeninger paper and our data
  intersect_data_upregulation <- list(
    "our_data_UpReg" = Gene_Names_MES_positive,
    "Groeningen_data" = unlist(ADR_MES_genes %>% filter(type == "MES") %>% dplyr::select(gene))
  )
  intersect_data_dwnregulation <- list(
    "our_data_DownReg" = Gene_Names_MES_negative,
    "Groeningen_data" = unlist(ADR_MES_genes %>% filter(type == "ADRN") %>% dplyr::select(gene))
  )

  p1 <- ggvenn::ggvenn(intersect_data_upregulation) + ggtitle(exp_name, " Overlap of CnR peaks (upreg. in MES)  with \n \t\t genes identified in Groeiningen paper, specific to MES \n")
  p2 <- ggvenn::ggvenn(intersect_data_dwnregulation) + ggtitle(exp_name, " Overlap of CnR peaks (upreg. in ADRN) with \n \t\t genes identified in Groeiningen paper, specific to ADRN \n")

  plot(p1)

  # perform GOs on the overlaps
  intersected_genes_upreg <- intersect(intersect_data_upregulation$Groeningen_data, intersect_data_upregulation$our_data_UpReg)
  background_genes_upreg <- length(unique(c(intersect_data_upregulation$Groeningen_data, intersect_data_upregulation$our_data_UpReg)))

  if (length(intersected_genes_upreg) != 0) {
    C5_GOBP_up <- hypeR::hypeR(
      signature = intersected_genes_upreg,
      genesets = gs_C5_GOBP,
      test = "hypergeometric",
      background = background_genes_upreg
    )

    C5_GOCC_up <- hypeR::hypeR(
      signature = intersected_genes_upreg,
      genesets = gs_C5_GOCC,
      test = "hypergeometric",
      background = background_genes_upreg
    )

    C5_GOMF_up <- hypeR::hypeR(
      signature = intersected_genes_upreg,
      genesets = gs_C5_GOMF,
      test = "hypergeometric",
      background = background_genes_upreg
    )

    C2_kegg_up <- hypeR::hypeR(
      signature = intersected_genes_upreg,
      genesets = gs_C2_kegg,
      test = "hypergeometric",
      background = background_genes_upreg
    )

    C2_reactome_up <- hypeR::hypeR(
      signature = intersected_genes_upreg,
      genesets = gs_C2_reactome,
      test = "hypergeometric",
      background = background_genes_upreg
    )

    plot(hypeR::hyp_dots(C5_GOBP_up, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOBP: upreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C5_GOCC_up, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOCC: upreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C5_GOMF_up, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOMF: upreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C2_kegg_up, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "KEGG: upreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C2_reactome_up, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "REACTOME: upreg in MES") + theme_bw())
  }

  plot(p2)

  intersected_genes_dwnreg <- intersect(intersect_data_dwnregulation$Groeningen_data, intersect_data_dwnregulation$our_data_DownReg)
  background_genes_dwnreg <- length(unique(c(intersect_data_dwnregulation$Groeningen_data, intersect_data_dwnregulation$our_data_DownReg)))

  if (length(intersected_genes_dwnreg) != 0) {
    C5_GOBP_dwn <- hypeR::hypeR(
      signature = intersected_genes_dwnreg,
      genesets = gs_C5_GOBP,
      test = "hypergeometric",
      background = background_genes_dwnreg
    )
    C5_GOCC_dwn <- hypeR::hypeR(
      signature = intersected_genes_dwnreg,
      genesets = gs_C5_GOCC,
      test = "hypergeometric",
      background = background_genes_dwnreg
    )
    C5_GOMF_dwn <- hypeR::hypeR(
      signature = intersected_genes_dwnreg,
      genesets = gs_C5_GOMF,
      test = "hypergeometric",
      background = background_genes_dwnreg
    )
    C2_kegg_dwn <- hypeR::hypeR(
      signature = intersected_genes_dwnreg,
      genesets = gs_C2_kegg,
      test = "hypergeometric",
      background = background_genes_dwnreg
    )
    C2_reactome_dwn <- hypeR::hypeR(
      signature = intersected_genes_dwnreg,
      genesets = gs_C2_reactome,
      test = "hypergeometric",
      background = background_genes_dwnreg
    )
    plot(hypeR::hyp_dots(C5_GOBP_dwn, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOBP: dwnreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C5_GOCC_dwn, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOCC: dwnreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C5_GOMF_dwn, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOMF: dwnreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C2_kegg_dwn, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "KEGG: dwnreg in MES") + theme_bw())
    plot(hypeR::hyp_dots(C2_reactome_dwn, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "REACTOME: dwnreg in MES") + theme_bw())
  }

  # check the total Gos for the CnR Upreg. data

  unique(intersect_data_upregulation$our_data_UpReg)

  C5_GOBP_up_total <- hypeR::hypeR(
    signature = unique(intersect_data_upregulation$our_data_UpReg),
    genesets = gs_C5_GOBP,
    test = "hypergeometric",
    background = background_genes_upreg
  )

  C5_GOCC_up_total <- hypeR::hypeR(
    signature = unique(intersect_data_upregulation$our_data_UpReg),
    genesets = gs_C5_GOCC,
    test = "hypergeometric",
    background = background_genes_upreg
  )

  C5_GOMF_up_total <- hypeR::hypeR(
    signature = unique(intersect_data_upregulation$our_data_UpReg),
    genesets = gs_C5_GOMF,
    test = "hypergeometric",
    background = background_genes_upreg
  )

  C2_kegg_up_total <- hypeR::hypeR(
    signature = unique(intersect_data_upregulation$our_data_UpReg),
    genesets = gs_C2_kegg,
    test = "hypergeometric",
    background = background_genes_upreg
  )

  C2_reactome_up_total <- hypeR::hypeR(
    signature = unique(intersect_data_upregulation$our_data_UpReg),
    genesets = gs_C2_reactome,
    test = "hypergeometric",
    background = background_genes_upreg
  )

  plot(hypeR::hyp_dots(C5_GOBP_up_total, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOBP: upreg in MES total") + theme_bw())
  plot(hypeR::hyp_dots(C5_GOCC_up_total, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOCC: upreg in MES total ") + theme_bw())
  plot(hypeR::hyp_dots(C5_GOMF_up_total, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "GOMF: upreg in MES total ") + theme_bw())
  plot(hypeR::hyp_dots(C2_kegg_up_total, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "KEGG: upreg in MES total ") + theme_bw())
  plot(hypeR::hyp_dots(C2_reactome_up_total, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = "REACTOME: upreg in MES total") + theme_bw())

  dev.off()
}



# load necessary datasets
# read xlxs with RNA-seq data
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/20230703/cell_type_MES_vs_ADR.xlsx", sheet = 1)
# read xlsx with ATACseq data
ATAC_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/20230510/DAR_results_3kb_promoters_M_VS_A.xlsx", sheet = 1)
# Subset  <-  to YAP, TAZ and Jun

DBobj_list_TAZ_YAP_Jun <- DBobj_list[c("YAP", "TAZ", "Jun")]

for (dfbobj in DBobj_list_TAZ_YAP_Jun) {
  factor_NAME <- dfbobj$DESeq2$DEdata$Factor[1]

  # load the object and extract the gene names
  rep <- as.data.frame(dba.report(dfbobj, contrast = 2))
  selector_bool <- as.data.frame(mcols(dfbobj$DESeq2$DEdata))$row.names.mcols_data. %in% rownames(rep)
  YAP_affected_genes <- unique(as.data.frame(mcols(dfbobj$DESeq2$DEdata))[selector_bool, ]$Gene.Name)

  # the list with all 3 datasets
  OVERLAP_list <- list(
    CnR_data = YAP_affected_genes,
    RNA = RNA_SEQ_data$gene_symbol,
    ATAC = ATAC_SEQ_data$gencode_gene_name
  )

  # open pdf
  pdf(file = paste0("~/workspace/neuroblastoma/results/20230919/", factor_NAME, "_CnR_RNA_ATAC_overlap.pdf"))

  # plot Venn diagram
  plot(ggvenn::ggvenn(
    OVERLAP_list,
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
    stroke_size = 0.5, set_name_size = 4
  ))

  # take the overlapping genes from 3 data sets and try to see what is the GOs
  deg_results <- intersect(intersect(OVERLAP_list$CnR_data, OVERLAP_list$RNA), OVERLAP_list$ATAC)
  write.table(deg_results, file = paste0("~/workspace/neuroblastoma/results/20230919/", factor_NAME, "_CnR_RNA_ATAC_overlap.tsv"), quote = F, col.names = F, row.names = F, sep = "\t")

  C5_GOBP <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOBP,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C5_GOCC <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOCC,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C5_GOMF <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOMF,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C2_kegg <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C2_kegg,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C2_reactome <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C2_reactome,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  plot(hypeR::hyp_dots(C5_GOBP, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOBP: ", factor_NAME, "_CnR_RNA_ATAC_overlap")) + theme_bw())
  plot(hypeR::hyp_dots(C5_GOCC, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOCC: ", factor_NAME, "_CnR_RNA_ATAC_overlap")) + theme_bw())
  plot(hypeR::hyp_dots(C5_GOMF, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOMF: ", factor_NAME, "_CnR_RNA_ATAC_overlap")) + theme_bw())
  plot(hypeR::hyp_dots(C2_kegg, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("KEGG: ", factor_NAME, "_CnR_RNA_ATAC_overlap")) + theme_bw())
  plot(hypeR::hyp_dots(C2_reactome, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("RACTOME: ", factor_NAME, "_CnR_RNA_ATAC_overlap")) + theme_bw())

  # compare this set to Groningen data - only 51 gene overlaps with Groningen data
  deg_results <- intersect(deg_results, ADR_MES_genes$gene)
  write.table(deg_results, file = paste0("~/workspace/neuroblastoma/results/20230919/", factor_NAME, "CnR_RNA_ATAC_overlap_with_Groningen_data.tsv"), quote = F, col.names = F, row.names = F, sep = "\t")

  C5_GOBP <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOBP,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C5_GOCC <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOCC,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C5_GOMF <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C5_GOMF,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C2_kegg <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C2_kegg,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  C2_reactome <- hypeR::hypeR(
    signature = deg_results,
    genesets = gs_C2_reactome,
    test = "hypergeometric",
    background = length(unique(unlist(OVERLAP_list)))
  )
  plot(hypeR::hyp_dots(C5_GOBP, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOBP: ", factor_NAME, "_CnR_RNA_ATAC_overlap_with_Groningen_data")) + theme_bw())
  plot(hypeR::hyp_dots(C5_GOCC, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOCC: ", factor_NAME, "_CnR_RNA_ATAC_overlap_with_Groningen_data")) + theme_bw())
  plot(hypeR::hyp_dots(C5_GOMF, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("GOMF: ", factor_NAME, "_CnR_RNA_ATAC_overlap_with_Groningen_data")) + theme_bw())
  plot(hypeR::hyp_dots(C2_kegg, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("KEGG: ", factor_NAME, "_CnR_RNA_ATAC_overlap_with_Groningen_data")) + theme_bw())
  plot(hypeR::hyp_dots(C2_reactome, merge = TRUE, fdr = 0.05, top = 20, abrv = 70, val = "fdr", title = paste0("RACTOME: ", factor_NAME, "_CnR_RNA_ATAC_overlap_with_Groningen_data")) + theme_bw())

  dev.off()
}



###########
peakAnno <- annotatePeak(files[[4]],
                         tssRegion = c(-3000, 3000),
                         TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"
)















###############
# try to identify active enh, promoters etc
# 4 categories -
# Promoters - within ±3 kb from TSS, or by annotation
# Active enh - overlap with H3k3me1 & H3k27ac
# inactive enh - overlap with h3k27ac








# Identify TFBS
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
motifsToScan_TEAD <- sapply(motifsToScan@listData, function(x){ str_detect(string = x@"name", pattern = "TEAD.") })
motifsToScan_TEAD <- motifsToScan[motifsToScan_TEAD]





# Using different subsets of peaks to check 
subsets_to_run <- c("TAZ", "YAP")

test_df <- SummarizedExperiment(assays = assay(DBobj_list[["TAZ"]][["DESeq2"]][["DEdata"]], rowRanges = DBobj_list[["TAZ"]][["merged"]])    )                    


for (sbst in subsets_to_run){
  
  example_counts <- chromVAR::addGCBias(object = test_df,
                                        genome = BSgenome.Hsapiens.UCSC.hg38)
  
  
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
                                            subject = DBobj_list[["TAZ"]][["merged"]], 
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
  ggsave(chrom_access_variability_plot, filename = file.path(res_dir, paste0("ChromVAR_chrom_access_variability_corr_plot_", sbst, ".pdf")))
  
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
  
  color.scheme <- c("#001219","#005F73", "#0A9396", "#94D2BD", "#E9D8A6", "#EE9B00", "#CA6702", "#BB3E03","#AE2012")
  
  heatmap <- pheatmap::pheatmap(as.matrix(devToPlot),
                                color = color.scheme,
                                #color = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
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
  
  ggsave(filename = file.path(res_dir, paste0("heatmap_ADRN_vs_MES_Motifs", sbst, "- used", ".png")), 
         plot=heatmap,
         width = 18, height = 20, units = "cm")
  
  ggsave(filename = file.path(res_dir, paste0("heatmap_ADRN_vs_MES_Motifs", sbst, "- used", ".eps")), 
         plot=heatmap,
         width = 18, height = 20, units = "cm")
}











### Analysis part
# Calculate for every factor an Enrichment score - very simplified and not reliable probably.
# load a table with genes that are specific to ADRN and MES cell lines + load annotation
mes_adrn_gene_list <- read.table("~/workspace/neuroblastoma/resources/mes_adrn_genes.tsv")
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(
  attributes = c(
    "hgnc_symbol", "chromosome_name",
    "start_position", "end_position"
  ),
  filters = c("hgnc_symbol"),
  values = mes_adrn_gene_list$V1,
  mart = mart
)
# remove scaffolds, keep only conventional chromosomes
results <- results[sapply(results$chromosome_name, nchar) < 3, ]
# keep the positional info and what cell type a gene belongs to
results <- merge.data.frame(mes_adrn_gene_list,
                            results,
                            by.y = "hgnc_symbol",
                            by.x = "V1"
)
Gresults <- makeGRangesFromDataFrame(results,
                                     keep.extra.columns = T,
                                     seqnames.field = "chromosome_name",
                                     start.field = "start_position",
                                     end.field = "end_position"
)

# subset GRanges for ADRN and MES genes and calculate scaling factor (to normalize the ES score) - universal for all targets
Gresults_ADRN <- Gresults[Gresults$V2 == "ADRN"]
ADRN_scaling_factor <- sum(as.data.frame(Gresults_ADRN)$end - as.data.frame(Gresults_ADRN)$start)
Gresults_MES <- Gresults[Gresults$V2 == "MES"]
MES_scaling_factor <- sum(as.data.frame(Gresults_MES)$end - as.data.frame(Gresults_MES)$start)

# now iterate over all proteins (YAP, TAZ, H3K4me1, H3K27ac, Jun)
ES_list <- data.frame(
  protein = character(0),
  ADRN = numeric(0),
  MES = numeric(0)
)



































###### will try the approach. Will try to build the rank only from promoters (±5kb around TSS, averaged affinity )
# MOdified function
# currently the cycle doesn't work - it breaks if there are not enough genes in a term
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES")) %>%
  dplyr::select(c(gene_symbol, Term))


for (DB_name in names(DBobj_list)) {
  DB_obj <- DBobj_list[[DB_name]]

  # filter only peaks that are located in promoters
  mcols(DB_obj$DESeq2$DEdata)$is_promoter <-
    abs(mcols(DB_obj$DESeq2$DEdata)$Distance.to.TSS) <= 5000
  # remove NA
  mcols(DB_obj$DESeq2$DEdata)$is_promoter[is.na(mcols(DB_obj$DESeq2$DEdata)$is_promoter)] <-
    FALSE
  tmp <-
    mcols(DB_obj$DESeq2$DEdata)[mcols(DB_obj$DESeq2$DEdata)$is_promoter, ]

  # tmp_report <- as.data.frame(dba.report(DBobj_list$Jun, contrast = 2))
  # tmp_report$row_name <- rownames(tmp_report)

  # tmp <- merge.data.frame(x = tmp, y = tmp_report, by.x = "row.names.mcols_data.", by.y = "row_name")

  tmp_gene_summary <- data.frame()
  for (gnname in unique(tmp$Gene.Name)) {
    indx <- gnname == tmp$Gene.Name
    tmp_gene_summary <- rbind(tmp_gene_summary, c(
      gnname,
      sum(indx),
      round(
        mean(tmp$Condition_M_vs_A[indx]),
        digits = 6
      )
    ))
  }

  colnames(tmp_gene_summary) <-
    c("GeneName", "num_of_genes", "MeanValue")
  tmp_gene_summary <- tmp_gene_summary %>% arrange(desc(MeanValue))
  tmp_gene_summary <- tmp_gene_summary[, c(1, 3)]

  tmp_vector <- as.numeric(tmp_gene_summary$MeanValue)
  names(tmp_vector) <- tmp_gene_summary$GeneName
  tmp_vector <- sort(tmp_vector, decreasing = TRUE)

  attr(tmp_vector, "names")

  # mes_adrn_gene_list <- read.table("~/workspace/neuroblastoma/resources/mes_adrn_genes.tsv")
  # mes_adrn_gene_list <- cbind(mes_adrn_gene_list$V2, mes_adrn_gene_list$V1)

  test <-
    clusterProfiler::GSEA(
      geneList = tmp_vector,
      TERM2GENE = mes_adrn_gene_list[mes_adrn_gene_list[, 1] == "MES", ],
      pvalueCutoff = 1,
      nPermSimple = 10000
    )
  p <- gseaplot3(test, geneSetID = 1, pvalue_table = T)
  p
  test <-
    clusterProfiler::GSEA(
      geneList = tmp_vector,
      TERM2GENE = mes_adrn_gene_list[mes_adrn_gene_list[, 1] == "ADRN", ],
      pvalueCutoff = 1,
      nPermSimple = 10000
    )
  p <- gseaplot3(test, geneSetID = 1, pvalue_table = T)
  p
}

# probably this result doesn't make a lot of sence. Or
#############################
# alternative approach -
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))


mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = "https://www.ensembl.org"
)

genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = mes_adrn_gene_list$ensembl_id,
  mart = mart
)


mes_adrn_gene_list <-
  merge.data.frame(
    x = mes_adrn_gene_list,
    y = genes,
    by.x = "ensembl_id",
    by.y = "ensembl_gene_id"
  ) %>%
  dplyr::select(Term, entrezgene_id) %>%
  dplyr::rename(gs_id = Term, gene_id = entrezgene_id)



##### Will try chipenrich

# create geneset ENTREZ ID
as.data.frame(TAZ_Peaks)

TAZ_Peaks_bed <- data.frame(
  seqnames = seqnames(TAZ_Peaks),
  start = start(TAZ_Peaks) - 1,
  end = end(TAZ_Peaks),
  names = c(rep(".", length(TAZ_Peaks))),
  scores = c(rep(".", length(TAZ_Peaks))),
  strands = strand(TAZ_Peaks)
)
TAZ_Peaks_bed[, 1] <- paste0("chr", TAZ_Peaks_bed[, 1])
TAZ_Peaks_bed <- TAZ_Peaks_bed %>% filter(seqnames != "chrGL000195.1")

results <- chipenrich(
  peaks = TAZ_Peaks_bed,
  genome = "hg38",
  genesets = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv",
  locusdef = "nearest_tss",
  qc_plots = FALSE,
  out_name = NULL,
  n_cores = 1
)

df <- data.frame(
  seqnames = seqnames(gr),
  starts = start(gr) - 1,
  ends = end(gr),
  names = c(rep(".", length(gr))),
  scores = c(rep(".", length(gr))),
  strands = strand(gr)
)




tmp@listData[["Condition_M_vs_A"]]

# need to pickup only peaks that belong to promoters
abs(DBobj_list$Jun[["DESeq2"]][["DEdata"]]@rowRanges@elementMetadata@listData[["Distance.to.TSS"]]) <= 3000








library(enrichplot)
gseaplot2(results, geneSetID = 1, title = results$Description[1])



gene_id <- results[["peaks"]][["gene_symbol"]][match(as.numeric(unlist(str_split(results[["results"]][["Geneset.Peak.Genes"]][1], pattern = ", "))), results[["peaks"]][["gene_id"]])]

results[["peaks"]][["gene_symbol"]][match(results[["peaks_per_gene"]][["gene_id"]], results[["peaks"]][["gene_id"]])]




##### try cluster profiler
library(clusterProfiler)
enricher(gene = unique(results[["peaks"]][["nearest_tss_symbol"]]), TERM2GENE = term_to_gene, minGSSize = 1)@result

GSEA(unique(results[["peaks"]][["nearest_tss_symbol"]]), TERM2GENE = term_to_gene, minGSSize = 1)



genes <- letters[1:15]
gs_df <- data.frame(
  "gs_name" = c(rep("genesetX", 10), rep("genesetY", 25)),
  "entrez_gene" = c(letters[1:10], letters[2:26])
)
enricher(gene = genes, TERM2GENE = gs_df, minGSSize = 1)@result


#### will try other approach
annoData <- ChIPpeakAnno::toGRanges(data = TxDb.Hsapiens.UCSC.hg38.knownGene)
peakdata <- DB
macs.anno <- ChIPpeakAnno::annotatePeakInBatch(GRanges(peakdata), AnnotationData = annoData)

peakAnno <- annotatePeak(files[[4]],
  tssRegion = c(-3000, 3000),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db"
)

annotatePeakInBatch()

### Chipenrich didn't work so well - i need to get plots.
### Will try other package
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
promoter <- getPromoters(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(GRanges(JUN_MES), windows = promoter)
peakHeatmap(GRanges(JUN_MES), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, upstream = 3000, downstream = 3000)

peakHeatmap(
  peak = GRanges(JUN_MES),
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  upstream = 3000,
  downstream = 3000
)


plotPeakProf2(
  peak = GRanges(JUN_MES), upstream = rel(0.2), downstream = rel(0.2),
  conf = 0.95, by = "gene", type = "body", nbin = 800,
  TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, ignore_strand = F
)

























# Will try to see the overlap between peaks of ATAC and CnR
# load ATAC_dds
ATAC_dds <- readRDS("~/workspace/neuroblastoma/RDSs/ATAC_dds.RDS")
ATAC_dds_CLB <- ATAC_dds[, ATAC_dds$cell_line == "CLBM"]

ATAC_dds_CLB@rowRanges@seqnames <- gsub(ATAC_dds_CLB@rowRanges@seqnames,
  pattern = "chr",
  replacement = ""
)
ATAC_ranges <- ATAC_dds_CLB@rowRanges
seqlevelsStyle(ATAC_ranges) <- "UCSC"

# load H3k27Ac
H3K27ac_CLBM <- dba.mask(DBobj_list[["H3K27ac"]],
  attribute = DBA_TISSUE,
  value = "CLB-Ma",
  bApply = T
)
H3K27ac_CLBM <- dba.contrast(H3K27ac_CLBM,
  minMembers = 2,
  reorderMeta = list(Condition = "A")
)
H3K27ac_CLBM <- dba.analyze(H3K27ac_CLBM)
H3K27ac_CLBM_ranges <- dba.peakset(H3K27ac_CLBM, bRetrieve = T)
seqlevelsStyle(H3K27ac_CLBM_ranges) <- "UCSC"

# load h3k4me1
H3K4me1_CLBM <- dba.mask(DBobj_list[["H3K4me1"]],
  attribute = DBA_TISSUE,
  value = "CLB-Ma",
  bApply = T
)
H3K4me1_CLBM <- dba.contrast(H3K4me1_CLBM,
  minMembers = 2,
  reorderMeta = list(Condition = "A")
)
H3K4me1_CLBM <- dba.analyze(H3K4me1_CLBM)
H3K4me1_CLBM_ranges <- dba.peakset(H3K4me1_CLBM, bRetrieve = T)
seqlevelsStyle(H3K4me1_CLBM_ranges) <- "UCSC"

ol <- ChIPpeakAnno::findOverlapsOfPeaks(ATAC_ranges, H3K27ac_CLBM_ranges, H3K4me1_CLBM_ranges)
ChIPpeakAnno::makeVennDiagram(ol,
  fill = c("#009E73", "#F0E442", "blue"), # circle fill color
  col = c("#D55E00", "#0072B2", "black"), # circle border color
  cat.col = c("#D55E00", "#0072B2", "black")
)

# load yap
YAP_CLBM <- dba.mask(DBobj_list[["YAP"]],
  attribute = DBA_TISSUE,
  value = "CLB-Ma",
  bApply = T
)
# YAP_CLBM <- dba.contrast(YAP_CLBM,
#                         minMembers=2,
#                         reorderMeta=list(Condition="A"))
# YAP_CLBM <- dba.analyze(YAP_CLBM)
YAP_CLBM_ranges <- dba.peakset(YAP_CLBM, bRetrieve = T)
seqlevelsStyle(YAP_CLBM_ranges) <- "UCSC"

# load TAZ
TAZ_CLBM <- dba.mask(DBobj_list[["TAZ"]],
  attribute = DBA_TISSUE,
  value = "CLB-Ma",
  bApply = T
)
# TAZ_CLBM <- dba.contrast(TAZ_CLBM,
#                         minMembers=2,
#                         reorderMeta=list(Condition="A"))
# TAZ_CLBM <- dba.analyze(TAZ_CLBM)
TAZ_CLBM_ranges <- dba.peakset(TAZ_CLBM, bRetrieve = T)
seqlevelsStyle(TAZ_CLBM_ranges) <- "UCSC"

ol <- ChIPpeakAnno::findOverlapsOfPeaks(YAP_CLBM_ranges, H3K27ac_CLBM_ranges, TAZ_CLBM_ranges)

ChIPpeakAnno::makeVennDiagram(ol,
  fill = c("#009E73", "#F0E442", "blue"), # circle fill color
  col = c("#D55E00", "#0072B2", "black"), # circle border color
  cat.col = c("#D55E00", "#0072B2", "black")
)

# 974 peaks are overlapping between all 3 sets. Try to annotate
YAP_TAZ_H3K27ac_Peaks <- ol$peaklist$`YAP_CLBM_ranges///H3K27ac_CLBM_ranges///TAZ_CLBM_ranges`
write.table(cbind(as.data.frame(YAP_TAZ_H3K27ac_Peaks)[, 1:3], "."),
  file = "~/workspace/neuroblastoma/tmp/YAP_TAZ_H3k27ac_Peaks_overlap.bed",
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = F
)

YAP_TAZ_H3K27ac_Peaks_annot <- read.csv2("~/workspace/neuroblastoma/tmp/YAP_TAZ_H3k27ac_Peaks_overlap_annotated.txt", sep = "\t")
# deg_results <- YAP_TAZ_H3K27ac_Peaks_annot$Gene.Name


annoData <- ChIPpeakAnno::toGRanges(data = TxDb.Hsapiens.UCSC.hg38.knownGene)
macs.anno <- annotatePeakInBatch(YAP_TAZ_ATAC_Peaks, AnnotationData = annoData)

library(org.Hs.eg.db)
macs.anno <- addGeneIDs(
  annotatedPeak = macs.anno,
  orgAnn = "org.Hs.eg.db",
  IDs2Add = "symbol"
)

deg_results <- unique(macs.anno@elementMetadata@listData[["symbol"]])


overlaps.anno$gene_name <- annoData$gene_name[match(
  overlaps.anno$feature,
  names(annoData)
)]
head(overlaps.anno)




subsetByOverlaps(ATAC_ranges, H3K27ac_CLBM_ranges)

library(ChIPpeakAnno)

ol <- findOverlapsOfPeaks(YAP_CLBM_ranges, H3K27ac_CLBM_ranges, TAZ_CLBM_ranges)

ol2 <- findOverlapsOfPeaks(H3K27ac_CLBM_ranges, H3K4me1_CLBM_ranges)

makeVennDiagram(ol,
  fill = c("#009E73", "#F0E442", "blue"), # circle fill color
  col = c("#D55E00", "#0072B2", "black"), # circle border color
  cat.col = c("#D55E00", "#0072B2", "black")
)

makeVennDiagram(ol2)

countOverlaps(dba.peakset(H3K27ac_CLBM, bRetrieve = T),
  ATAC_dds_CLB@rowRanges,
  minoverlap = 1
)

as.data.frame(dba.peakset(H3K27ac_CLBM, bRetrieve = T)) %>%
  dplyr::select(columns = c("seqnames", "start", "end")) %>%
  write.table(., file = "./neuroblastoma/results/H3K27ac_export.bed", sep = "\t", quote = F, col.names = F, row.names = F)

as.data.frame(ATAC_dds_CLB@rowRanges) %>%
  dplyr::select(columns = c("seqnames", "start", "end")) %>%
  write.table(., file = "./neuroblastoma/results/atac_export.bed", sep = "\t", quote = F, col.names = F, row.names = F)





# rep <- rep[abs(rep$Fold) > 2,]
#profiles <- dba.plotProfile(dfbobj,
#                             merge = c(DBA_TISSUE, DBA_REPLICATE),
#                            contrast=2,
#                           sites = repList)
#dba.plotProfile(profiles )

# profiles <- dba.plotProfile(dfbobj, merge = c(DBA_REPLICATE), contrast=1)
# dba.plotProfile(profiles )




rep <- dba.report(dfbobj, contrast = 2)

repUP <- rep[rep$Fold > 2, ]
repDWN <- rep[rep$Fold < -1, ]

repUP <- repUP[order(repUP$Fold, decreasing = TRUE), ]
repDWN <- repDWN[order(repDWN$Fold, decreasing = FALSE), ]

repList <- GRangesList(
  UP = repUP,
  DWN = repDWN
)

rep <- rep[abs(rep$Fold) > 2, ]

profiles <- dba.plotProfile(dfbobj,
  merge = c(DBA_TISSUE, DBA_REPLICATE),
  contrast = 2,
  sites = repList
)
dba.plotProfile(profiles)

sum(dfbobj.DB$Fold < 0)


tmp <- dba.plotProfile(dfbobj, sites = repList, contrast = 2)





# dba.peakset(dfbobj, -DBA_REPLICATE)

dfbobj <- dba.normalize(dfbobj)

dfbobj <- dba.contrast(dfbobj,
  minMembers = 2,
  reorderMeta = list(Condition = "A")
)

dfbobj <- dba.analyze(dfbobj)


dba.show(dfbobj, bContrasts = TRUE)
plot(dfbobj, contrast = 1)
plot(dfbobj, contrast = 2)

dfbobj.DB <- dba.report(dfbobj, contrast = 2)






dba.plotPCA(dfbobj, contrast = 2)

dba.plotMA(dfbobj)

dba.plotVolcano(dfbobj, contrast = 2)


hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
readscores <- dba.plotHeatmap(dfbobjsubsample,
  contrast = 2, correlations = FALSE,
  scale = "row", colScheme = hmap, bUsePval = T
)






# load YAP_TAZ_enhancers
yap_taz_enh <- read.table("./neuroblastoma/resources/YAP_TAZ_affected_enhancers.tsv", sep = "\t", dec = ".", header = TRUE)
yap_taz_enh <- GRanges(yap_taz_enh)
yap_taz_enh

gr <- subsetByOverlaps(yap_taz_enh, GRanges(dfbobj$peaks[[5]]))
df <- data.frame(
  seqnames = seqnames(gr),
  starts = start(gr) - 1,
  ends = end(gr),
  names = c(gr$TARGET_GENE),
  scores = c(rep(".", length(gr))),
  strands = strand(gr)
)

write.table(df, file = "foo.bed", quote = F, sep = "\t", row.names = F, col.names = F)





findOverlaps(yap_taz_enh, GRanges(dfbobj$peaks[[5]]))
subsetByOverlaps(yap_taz_enh, GRanges(dfbobj$peaks[[5]]))

count_lines(files)



####################



dfbobjsubsampleADRN <- dba.count(dfbobj, peaks = subset_peak_set_ADRN)
dfbobjsubsample <- dba.normalize(dfbobjsubsample)
dfbobjsubsample <- dba.contrast(dfbobjsubsample,
  minMembers = 2,
  reorderMeta = list(Condition = "A")
)
dfbobjsubsample <- dba.analyze(dfbobjsubsample)







subset_peak_set <- subsetByOverlaps(Gresults, subset_peak_set)

dfbobjsubsample <- dba.count(dfbobj, peaks = subset_peak_set)


dfbobjsubsample <- dba.normalize(dfbobjsubsample)
dfbobjsubsample <- dba.contrast(dfbobjsubsample,
  minMembers = 2,
  reorderMeta = list(Condition = "A")
)
dfbobjsubsample <- dba.analyze(dfbobjsubsample)




subsetByOverlaps(dfbobj, subset_peak_set)












YAP_conssensus <- as.data.frame(dfbobj[["merged"]])
YAP_conssensus <- cbind(
  paste0(YAP_conssensus$CHR, "_", YAP_conssensus$START, "_", YAP_conssensus$END),
  YAP_conssensus,
  "."
)
YAP_conssensus <- cbind(
  row.names(YAP_conssensus),
  YAP_conssensus,
  "."
)


YAP_conssensus$CHR <- paste0("chr", YAP_conssensus$CHR)

write.table(YAP_conssensus,
  file = "neuroblastoma/results/20230828/consensus_peaks/YAP.bed",
  sep = "\t",
  row.names = F,
  col.names = F,
  quote = F
)

mcols_data <- mcols(dfbobj$DESeq2$DEdata)
mcols_data <- cbind(row.names(mcols_data), mcols_data)

YAP_annot <- read.csv2("./neuroblastoma/results/20230828/consensus_peaks/yap_output.txt", sep = "\t")

mcols(dfbobj$DESeq2$DEdata) <- merge.data.frame(mcols_data, YAP_annot, by.x = 1, by.y = 1, sort = F)
