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
if(!dir.exists("/home/rstudio/.cache/R/ExperimentHub")) {dir.create("/home/rstudio/.cache/R/ExperimentHub")}

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
  "chipEnrichAndExport",
  "createTermsTable",
  "subsetPeaksInSamples",
  "LolipopEnrichmentPlot"
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

options(timeout = 30000)
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

# this is necessary to prevent VennDiagrap package from spamming logs into workspace folder
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

system("unzip neuroblastoma/data/CnR/peaks.zip -d neuroblastoma/data/CnR/peaks/")

system("zip -F neuroblastoma/data/CnR/consensus_peaks.zip --out neuroblastoma/data/CnR/consensus_peaks_u.zip")
system("unzip neuroblastoma/data/CnR/consensus_peaks_u.zip -d neuroblastoma/data/CnR/")
system("rm neuroblastoma/data/CnR/consensus_peaks_u.zip")

data_folder <- "/home/rstudio/workspace/neuroblastoma/data/CnR"

files <- unlist(list.files(path = file.path(data_folder, "consensus_peaks"), pattern = "dba", full.names = TRUE))
sample_names <- stringr::str_extract(string = files, 
                                     pattern = ".*(dba_)(.*).RData", group = 2)
DBobj_list <- list()
# DBobj_list_cons <- list()
for (target_prot in unique(sample_names)) {
  print(paste("File", target_prot, "exists, loading ..."))
  fl_path <- file.path(data_folder, paste0("consensus_peaks/dba_", target_prot))
  dfbobj <- dba.load(fl_path, dir = "", pre = "")
  DBobj_list[[target_prot]] <- dfbobj
#  DBobj_list_cons[[target_prot]] <- dba.peakset(dfbobj, bRetrieve = TRUE)
}



# ######
# dba.overlap(DBobj_list$H3K27ac, mode = DBA_OLAP_RATE)
# dba.overlap(DBobj_list$H3K4me1, mode = DBA_OLAP_RATE)
# dba.overlap(DBobj_list$Jun, mode = DBA_OLAP_RATE)
# #####

# # Load all peaks from .narrowPeaks files
folder_path <- "~/workspace/neuroblastoma/data/CnR/peaks/"
files <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.narrowPeak$", full.names = TRUE))
summits <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.macs2_summits.bed$", full.names = TRUE))

# assemble the samples object
samples <- data.frame(list(SampleID = str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R(\\d)")))
samples$Tissue <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)")
samples$Factor <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R", group = 3)
samples$Condition <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)", group = 2)
samples$Replicate <- str_extract(files, pattern = "_R(\\d)", group = 1)
samples$Peaks <- files
samples$PeakCaller <- "narrow"
samples$PeakSummits <- summits
peaks_count <- c()
peaks_count <- sapply(files, function(x) genio::count_lines(x) - 1)
samples$Peaks_count <- peaks_count

# Loading all narrowPeaks
DBobj_list_cons <- list()
for (target_prot in unique(samples$Factor)) {
  print(c("Current protein is ", target_prot))
  dfbobj <- dba(sampleSheet = samples %>% filter(., Factor == target_prot))
  DBobj_list_cons[[target_prot]] <- dba.peakset(dfbobj, bRetrieve = TRUE)
}

peakAnnoList <- list()
# check features distribution
for (factor in names(DBobj_list_cons)){
  test_region <- DBobj_list_cons[[factor]]
  seqlevelsStyle(test_region) <- "UCSC"
  peakAnnoList[[factor]] <- ChIPseeker::annotatePeak(peak = test_region, 
                                         tssRegion = c(-3000, 3000),
                                         TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                         annoDb = "org.Hs.eg.db", 
                                         level = "gene")
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
mes_adrn_gene_list <- createTermsTable(path_to_RNAseq_xlsx_table = "~/workspace/neuroblastoma/results/RNA-seq/cell_type_MES_vs_ADR.xlsx", 
                                       name_for_negative_values = "ADRN",
                                       name_for_positive_values = "MES",
                                       mart = mart
)
write.table(
  x = mes_adrn_gene_list,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)

# Create a term from K975 experiment 24h
RNA_K975_SEQ_data_24h <- createTermsTable(path_to_RNAseq_xlsx_table = "~/workspace/neuroblastoma/results/RNA-seq_yap_taz_inhibition/cell_type_24h_vs_control.xlsx",
                                          name_for_negative_values = "K975_24h_down",
                                          name_for_positive_values = "K975_24h_up",
                                          mart = mart
)
write.table(
  x = RNA_K975_SEQ_data_24h,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_24h.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)


# Create a term from K975 experiment 48h
RNA_K975_SEQ_data_48h <- createTermsTable(path_to_RNAseq_xlsx_table = "~/workspace/neuroblastoma/results/RNA-seq_yap_taz_inhibition/cell_type_48h_vs_control.xlsx",
                                          name_for_negative_values = "K975_48h_down",
                                          name_for_positive_values = "K975_48h_up",
                                          mart = mart
)
write.table(
  x = RNA_K975_SEQ_data_48h,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_48h.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)


RNA_K975_SEQ_data_24h_48h <- rbind(
  data.frame(gs_id = "K975_24h_48h_up",
             gene_id = intersect(RNA_K975_SEQ_data_24h %>% filter(gs_id == "K975_24h_up") %>% pull(gene_id),
                                 RNA_K975_SEQ_data_48h %>% filter(gs_id == "K975_48h_up") %>% pull(gene_id))
  ),
  data.frame(gs_id = "K975_24h_48h_down",
             gene_id = intersect(RNA_K975_SEQ_data_24h %>% filter(gs_id == "K975_24h_down") %>% pull(gene_id),
                                 RNA_K975_SEQ_data_48h %>% filter(gs_id == "K975_48h_down") %>% pull(gene_id))
  )
)
write.table(
  x = RNA_K975_SEQ_data_24h_48h,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_24h_48h.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)

# Create a term from JunDN experiment 
RNA_JunDN_SEQ_data <- createTermsTable(path_to_RNAseq_xlsx_table = "~/workspace/neuroblastoma/results/RNA-seq_yap_taz_inhibition/Jun_DN/JunDN_dox_vs_ctrl.xlsx",
                                          name_for_negative_values = "JunDN_down",
                                          name_for_positive_values = "JunDN_up",
                                          mart = mart
)
write.table(
  x = RNA_JunDN_SEQ_data,
  file = "~/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_JunDN_SEQ_data.tsv",
  sep = "\t",
  row.names = F,
  quote = F
)

locusdef <-  "5kb"

#our_rna_seq_terms <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv"
#our_rna_seq_terms <- mes_adrn_gene_list

res_dir <- "/home/rstudio/workspace/neuroblastoma/results/CnR/"
res_dir_k975 <- "/home/rstudio/workspace/neuroblastoma/results/CnR/k975/"

#our_rna_seq_terms <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv"
#our_rna_seq_terms <- mes_adrn_gene_list
our_rna_seq_terms <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_seq.tsv"
our_rna_seq_terms_k975_24h <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_24h.tsv"
our_rna_seq_terms_k975_48h <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_48h.tsv"
our_rna_seq_terms_k975_24h_48h <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_k975_RNA_seq_24h_48h.tsv"
our_rna_seq_terms_JunDN <- "/home/rstudio/workspace/neuroblastoma/resources/mes_adrn_GS_frm_RNA_JunDN_SEQ_data.tsv"

rna_seq_terms_list <- list(
  rna_seq_terms = list(name = "Our_data_", path = our_rna_seq_terms),
  k975_24h = list(name = "Our_k975_24h_data_", path = our_rna_seq_terms_k975_24h),
  k975_48h = list(name = "Our_k975_48h_data_", path =our_rna_seq_terms_k975_48h),
  k975_24h_48h = list(name = "Our_k975_24h_48h_data_", path = our_rna_seq_terms_k975_24h_48h),
  JunDN = list(name = "Our_JunDN_data_", path = our_rna_seq_terms_JunDN)
)


for(TF_name in names(DBobj_list)){
  #TF_name <- "YAP"  
  print(paste0("Processing ", TF_name))
  
  # subset peaks that up-regulated in MES samples
  p_MES_up <- subsetPeaksInSamples(DBobj_list, TF_name, select_up_or_down_regulated = "up", remove_scaffolds = TRUE)
  p_MES_dwn <- subsetPeaksInSamples(DBobj_list, TF_name, select_up_or_down_regulated = "down", remove_scaffolds = TRUE)
  #plot_dist_to_tss(peaks = p_MES_dwn, genome = "hg38")
  
  if (nrow(p_MES_up) > 0) {
    for(rna_seq_term in rna_seq_terms_list){
      chipEnrichAndExport(peaks = p_MES_up,
                          peaksName = "MES", 
                          TF_name = TF_name, 
                          res_dir = res_dir, 
                          genesets = rna_seq_term$path, 
                          genesets_name = rna_seq_term$name,
                          locusdef = locusdef
      )
    }
  }
  
  if (nrow(p_MES_dwn) > 0){
    for(rna_seq_term in rna_seq_terms_list){
      chipEnrichAndExport(peaks = p_MES_dwn,
                          peaksName = "ADRN", 
                          TF_name = TF_name, 
                          res_dir = res_dir, 
                          genesets = rna_seq_term$path, 
                          genesets_name = rna_seq_term$name,
                          locusdef = locusdef
      )
    }
  }
}

# Addition for a reviewer - combining YAP/TAZ and YAP/TAZ/Jun
YAP_p_MES_up <- subsetPeaksInSamples(DBobj_list, "YAP", contrast = 2, select_up_or_down_regulated = "up", remove_scaffolds = TRUE)
TAZ_p_MES_up <- subsetPeaksInSamples(DBobj_list, "TAZ", contrast = 2, select_up_or_down_regulated = "up", remove_scaffolds = TRUE)
JUN_p_MES_up <- subsetPeaksInSamples(DBobj_list, "Jun", contrast = 2, select_up_or_down_regulated = "up", remove_scaffolds = TRUE)
# process up-regulated peaks
# remove scaffolds - their name is longer than 2 charcaters
ol <- ChIPpeakAnno::findOverlapsOfPeaks(GRanges(YAP_p_MES_up), GRanges(TAZ_p_MES_up))
YAP_TAZ_MES_p_up <- as.data.frame(ol$mergedPeaks)
ol <- ChIPpeakAnno::findOverlapsOfPeaks(GRanges(YAP_TAZ_MES_p_up), GRanges(JUN_p_MES_up))
YAP_TAZ_JUN_MES_p_up <- as.data.frame(ol$mergedPeaks)

for(rna_seq_term in rna_seq_terms_list){
  chipEnrichAndExport(peaks = YAP_TAZ_MES_p_up,
                      peaksName = "MES", 
                      TF_name = "yap-taz", 
                      res_dir = res_dir, 
                      genesets = rna_seq_term$path, 
                      genesets_name = rna_seq_term$name,
                      locusdef = locusdef
  )
}

for(rna_seq_term in rna_seq_terms_list){
  chipEnrichAndExport(peaks = YAP_TAZ_JUN_MES_p_up,
                      peaksName = "MES", 
                      TF_name = "yap-taz-jun", 
                      res_dir = res_dir, 
                      genesets = rna_seq_term$path, 
                      genesets_name = rna_seq_term$name,
                      locusdef = locusdef
  )
}


YAP_p_MES_down <- subsetPeaksInSamples(DBobj_list, "YAP", contrast = 2, select_up_or_down_regulated = "down", remove_scaffolds = TRUE)
TAZ_p_MES_down <- subsetPeaksInSamples(DBobj_list, "TAZ", contrast = 2, select_up_or_down_regulated = "down", remove_scaffolds = TRUE)
JUN_p_MES_down <- subsetPeaksInSamples(DBobj_list, "Jun", contrast = 2, select_up_or_down_regulated = "down", remove_scaffolds = TRUE)

ol <- ChIPpeakAnno::findOverlapsOfPeaks(GRanges(YAP_p_MES_down), GRanges(TAZ_p_MES_down))
YAP_TAZ_MES_p_down <- as.data.frame(ol$mergedPeaks)
ol <- ChIPpeakAnno::findOverlapsOfPeaks(YAP_TAZ_MES_p_down, GRanges(JUN_p_MES_down))
YAP_TAZ_JUN_MES_p_down <- as.data.frame(ol$mergedPeaks)

for(rna_seq_term in rna_seq_terms_list){
  chipEnrichAndExport(peaks = YAP_TAZ_MES_p_down,
                      peaksName = "ADRN", 
                      TF_name = "yap-taz", 
                      res_dir = res_dir, 
                      genesets = rna_seq_term$path, 
                      genesets_name = rna_seq_term$name,
                      locusdef = locusdef
  )
}

for(rna_seq_term in rna_seq_terms_list){
  chipEnrichAndExport(peaks = YAP_TAZ_JUN_MES_p_down,
                      peaksName = "ADRN", 
                      TF_name = "yap-taz-jun", 
                      res_dir = res_dir, 
                      genesets = rna_seq_term$path, 
                      genesets_name = rna_seq_term$name,
                      locusdef = locusdef
  )
}
  
####### visualize our data #########
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
pattern <-  ".*(H3K27ac|H3K4me1|Jun|TAZ|YAP|yap-taz|yap-taz-jun)_(Our_data|Our_k975_24h_data|Our_k975_48h_data|Our_k975_24h_48h_data|Our_JunDN_data).*(ADRN|MES).*"

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

pdf(file = paste0(res_dir, "enrichment_lolipop_plots.pdf"))
# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_data", Selected_peaks == "MES")
# making lolipop plot with enrichments
plot <- LolipopEnrichmentPlot(summary_table = summary_table_MES,
       title = "Enrichment of peaks in MES samples in MES/ADR-specific regions determined from RNA-seq",
       p.val_treshold = 0.05)
plot

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_data", Selected_peaks == "ADRN")
# making lolipop plot with enrichments
plot <- LolipopEnrichmentPlot(summary_table = summary_table_ADRN,
                              title = "Enrichment of peaks in ADRN samples in MES/ADR-specific regions determined from RNA-seq",
                              p.val_treshold = 0.05)
plot


# Visualization, but with genes that comes from K975 24 h inhibition assay
# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_k975_24h_data", Selected_peaks == "MES")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_MES,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in MES samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 24h"
                              )
plot

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_k975_24h_data", Selected_peaks == "ADRN")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_ADRN,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in ADRN samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 24h"
)
plot

# Visualization, but with genes that comes from 48h K975 inhibition assay
# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_k975_48h_data", Selected_peaks == "MES")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_MES,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in MES samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 48h"
)
plot

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_k975_48h_data", Selected_peaks == "ADRN")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_ADRN,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in ADRN samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 48h"
)
plot

# Visualization, but with genes that comes from 48h K975 inhibition assay
# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_k975_24h_48h_data", Selected_peaks == "MES")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_MES,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in MES samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 24_48h"
)
plot

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_k975_24h_48h_data", Selected_peaks == "ADRN")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_ADRN,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in ADRN samples in MES/ADR-specific \n regions determined from K-975 RNA-seq 24_48h"
)
plot

# Visualization, but with genes that comes from JunDN inhibition assay
# Make a plot for MES data
summary_table_MES <- summary_table %>% filter(Data_source == "Our_JunDN_data", Selected_peaks == "MES")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_MES,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in MES samples in MES/ADR-specific \n regions determined from Our JunDN data"
)
plot

# Make a plot for ADRN data
summary_table_ADRN <- summary_table %>% filter(Data_source == "Our_JunDN_data", Selected_peaks == "ADRN")
plot <- LolipopEnrichmentPlot(summary_table = summary_table_ADRN,
                              p.val_treshold = 0.05,
                              title = "Enrichment of peaks in ADRN samples in MES/ADR-specific \n regions determined from Our JunDN data"
)
plot


dev.off()

### Check distribution of features 
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

peakAnnoList <- list()
DBobj_list_to_annotate <- list(
  TAZ = TAZ_Peaks, 
  YAP = YAP_Peaks, 
  Jun = Jun_Peaks,
  YAP_TAZ = YAP_TAZ_peaks,
  YAP_TAZ_JUN = YAP_TAZ_JUN_peaks)
# check features distribution
for (factor in names(DBobj_list_to_annotate)){
  test_region <- DBobj_list_to_annotate[[factor]]
  seqlevelsStyle(test_region) <- "UCSC"
  peakAnnoList[[factor]] <- ChIPseeker::annotatePeak(peak = test_region, 
                                                     tssRegion = c(-3000, 3000),
                                                     TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                                     annoDb = "org.Hs.eg.db", 
                                                     level = "gene")
}
plotAnnoBar(peakAnnoList)

# Check what is the percentage on promoters, enhancers is in the peaks
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
  peaks <- reduce(YAP_TAZ_JUN_peaks_list[[peak_names]])
  
  active_promoters_percent <- length(findOverlaps(peaks, active_promoters))
  poised_promoters_percent <- length(findOverlaps(peaks, poised_promoters))
  active_enhancer_percent <- length(findOverlaps(peaks, active_enhancer))
  poised_enhancer_percent <- length(findOverlaps(peaks, poised_enhancer))
  
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
  
  
# Trying to calculate the observed/expected ratio for active enhancers, promoters, bivalent promoters, poised enhancers
genome_used <- getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
seqnames(genome_used) <- sub("^chr", "", seqnames(genome_used))
seqnames(genome_used) <- sub("^.{1,2}_", "", seqnames(genome_used))
seqnames(genome_used) <- sub("v", ".", seqnames(genome_used))
seqnames(genome_used) <- sub("^M", "MT", seqnames(genome_used))

# load mask regions 
blacklist_regions <- read.delim("~/workspace/neuroblastoma/resources/hg38-blacklist.v2.bed", header = FALSE)
blacklist_regions$V1 <- sub("^chr", "", blacklist_regions$V1)
blacklist_regions <- makeGRangesFromDataFrame(blacklist_regions, 
                                              seqnames.field = "V1", 
                                              start.field = "V2", 
                                              end.field = "V3")







piechart_randomised <- tibble()
for(peak_names in names(YAP_TAZ_JUN_peaks_list)){
  
  #peak_names <- "Jun_Peaks"
  
  initial_peak_set <- reduce(YAP_TAZ_JUN_peaks_list[[peak_names]])
  random_Peaks <- initial_peak_set
  random_Peaks@seqnames <- droplevels(random_Peaks@seqnames)
  random_Peaks <- regioneR::randomizeRegions(random_Peaks,
                                             allow.overlaps = FALSE,
                                             genome = genome_used,
                                             per.chromosome = TRUE,
                                             mask = blacklist_regions)
  random_Peaks <- as.data.frame(random_Peaks, row.names = NULL, optional = FALSE)
  random_Peaks$seqnames <- droplevels(random_Peaks$seqnames)
  random_Peaks <- makeGRangesFromDataFrame(random_Peaks)
  
  active_promoters_percent <- length(findOverlaps(random_Peaks, active_promoters))
  poised_promoters_percent <- length(findOverlaps(random_Peaks, poised_promoters))
  active_enhancer_percent <- length(findOverlaps(random_Peaks, active_enhancer))
  poised_enhancer_percent <- length(findOverlaps(random_Peaks, poised_enhancer))

  unchar <- length(initial_peak_set) - 
    active_promoters_percent - 
    poised_promoters_percent - 
    active_enhancer_percent - 
    poised_enhancer_percent

  piechart_randomised <- rbind(piechart_randomised,
                    
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
  

piechart_randomised <- piechart_randomised %>% group_by(group) %>% mutate(percent = round(number_of_peaks / sum(number_of_peaks) * 100))
  
piechart_randomised <- piechart_randomised  %>% mutate(percent = round(number_of_peaks / sum(number_of_peaks) * 100))
  
p <- ggplot(piechart_randomised) +
  aes(x = group, y = percent, fill = Peaks) +
  geom_col() +
  scale_fill_brewer(palette = "Accent", direction = 1) +
  theme_minimal() +
  geom_text(aes(label = paste0(percent, "%")),
            position = position_stack(vjust = 0.5), size = 6
  )
plot(p)

# save for IGV
rtracklayer::export.bed(object = random_Peaks,  "~/workspace/neuroblastoma/temp_results/BEDs/JUN_peaks_for_IGV_randomized.bed")
rtracklayer::export.bed(object = YAP_TAZ_JUN_peaks_list$Jun_Peaks, "~/workspace/neuroblastoma/temp_results/BEDs/JUN_peaks_for_IGV.bed")
rtracklayer::export.bed(object = promoter, "~/workspace/neuroblastoma/temp_results/BEDs/Promoters_for_IGV.bed")

x = GenomicDistributions::calcChromBinsRef(makeGRangesFromDataFrame(as.data.frame(YAP_TAZ_JUN_peaks_list$Jun_Peaks) %>% mutate(seqnames = paste0("chr", seqnames))), "hg38")
GenomicDistributions::plotChromBins(x)
x = GenomicDistributions::calcChromBinsRef(makeGRangesFromDataFrame(as.data.frame(random_Peaks) %>% mutate(seqnames = paste0("chr", seqnames))), "hg38")
GenomicDistributions::plotChromBins(x)


# Check that namings are the same:
piechart_randomised$Peaks == piechart$Peaks

# Do fisher exact test for observed/expected peaks
fisher_result_df <- tibble()

for(peak_group_name in unique(piechart_randomised$group)){
  total_number_of_peaks <- length(reduce(YAP_TAZ_JUN_peaks_list[[peak_group_name]]))
  
  for(feature_group_name in unique(piechart_randomised$Peaks)){
    peak_subset_number <- piechart %>% 
      dplyr::filter(group == peak_group_name & Peaks == feature_group_name) %>%
      pull(number_of_peaks)
    randomized_peak_subset_number <- piechart_randomised %>% 
      dplyr::filter(group == peak_group_name & Peaks == feature_group_name) %>%
      pull(number_of_peaks)
  
    fisher_matrix <- matrix(c(peak_subset_number, total_number_of_peaks - peak_subset_number,
                              randomized_peak_subset_number, total_number_of_peaks - randomized_peak_subset_number),
                            nrow = 2)
    tmp <- fisher.test(fisher_matrix)
    
    fisher_result_df <- rbind(fisher_result_df,
                              tibble(group = peak_group_name,
                                     Peaks = feature_group_name,
                                     odds_ratio = tmp$estimate,
                                     p_value = tmp$p.value))
   
    }
}

fisher_result_df$odds_ratio_log <- log10(fisher_result_df$odds_ratio)
fisher_result_df$p_value_mod <- fisher_result_df$p_value
fisher_result_df$p_value_mod[fisher_result_df$p_value_mod < 1e-10] <- 1e-10
fisher_result_df$p_value_log <- -log10(fisher_result_df$p_value_mod)

# Plot
ggplot(fisher_result_df, aes(x = Peaks, y = odds_ratio_log, group = group)) +
  geom_segment(aes(xend = Peaks, y = 0, yend = odds_ratio_log), color = "gray") +
  geom_point(aes(size = p_value_log), color = "blue") +
  facet_wrap(~group, scales = "free_x") +
  labs(
    title = "Lollipop Plot",
    x = "Peaks",
    y = "-log10(Odds)",
    size = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
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

