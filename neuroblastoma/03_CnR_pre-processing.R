folder_path <- "~/workspace/neuroblastoma/data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/04_called_peaks/macs2"

# List all files with *.narrowPeak extension in the folder and summits
files <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.narrowPeak$", full.names = TRUE))
summits <- unlist(list.files(path = folder_path, pattern = "(CLB-Ma|SK-N-SH).*.macs2_summits.bed$", full.names = TRUE))

# assemble the samples object
samples <- data.frame(list(SampleID = str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R(\\d)")))
samples$Tissue <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)")
samples$Factor <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R", group = 3)
samples$Condition <- str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)", group = 2)
samples$Replicate <- str_extract(files, pattern = "_R(\\d)", group = 1)
samples$bamReads <- paste0(
  "/home/rstudio/workspace/neuroblastoma/data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/",
  str_extract(files, pattern = "(CLB-Ma|SK-N-SH)-(A|M)_(.*)_R(\\d)"),
  ".target.markdup.sorted.bam"
)
samples$controlID <- paste0(str_extract(files, pattern = "((CLB-Ma|SK-N-SH)-(A|M))", group = 1), "_IgG_R", samples$Replicate)
samples$bamControl <- paste0(
  "/home/rstudio/workspace/neuroblastoma/data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/",
  samples$controlID,
  ".target.markdup.sorted.bam"
)
samples$bamControl <- str_replace_all(samples$bamControl, 
                                      pattern = "CLB-Ma-M_IgG_R2.target.markdup.sorted.bam",
                                      replacement = "CLB-Ma-M_IgG_R1.target.markdup.sorted.bam")
samples$Peaks <- files
samples$PeakCaller <- "narrow"
samples$PeakSummits <- summits
peaks_count <- c()
peaks_count <- sapply(files, function(x) count_lines(x) - 1)
samples$Peaks_count <- peaks_count