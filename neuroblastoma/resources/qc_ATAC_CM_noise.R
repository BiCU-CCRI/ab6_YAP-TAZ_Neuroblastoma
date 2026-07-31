.libPaths(c("/home/abykov/R/4.5", .libPaths()))
suppressMessages({library(SummarizedExperiment); library(DESeq2); library(dplyr); library(stringr); library(matrixStats)})
se <- readRDS("~/neuroblastoma/data/ATACseq/consensus_peaks.mLb.clN.rds")
keep_s <- c("CM_control_ATAC_S165169_REP1","CM_24h_ATAC_S165170_REP1","CM_48h_ATAC_S165167_REP1",
"SH_control_ATAC_S165166_REP1","SH_24h_ATAC_S165164_REP1","SH_48h_ATAC_S165165_REP1",
"SH_DMSO_ctrl_rep1_REP1","SH_DMSO_ctrl_rep2_REP1","SH_24h_rep1_REP1","SH_24h_rep2_REP1",
"SH_48h_rep1_REP1","SH_48h_rep2_REP1","CM_DMSO_ctrl_rep1_REP1","CM_DMSO_ctrl_rep2_REP1",
"CM_24h_rep1_REP1","CM_24h_rep2_REP1","CM_48h_rep1_REP1","CM_48h_rep2_REP1")
dds <- se[, se$sample %in% keep_s]
dds$line  <- ifelse(str_detect(dds$sample,"^CM"),"CM","SH")
dds$time  <- case_when(str_detect(dds$sample,"24h")~"24h",str_detect(dds$sample,"48h")~"48h",TRUE~"control")
dds$batch <- ifelse(str_detect(dds$sample,"_ATAC_S\\d+"),"legacy","new")
design(dds) <- ~1
keep <- rowSums(fpm(dds,robust=TRUE)>1) >= 2
dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
ln <- log2(counts(dds, normalized=TRUE) + 1)          # normalized log2 counts
colnames(ln) <- paste(dds$line, dds$time, dds$batch, sep="_")
# average the two new replicates into one profile per condition, keep legacy as-is
cond <- colnames(ln)
avg <- sapply(unique(cond), function(cc) rowMeans(ln[, cond==cc, drop=FALSE]))

cat("=====================================================================\n")
cat("A) CROSS-BATCH AGREEMENT: does the NEW batch reproduce the LEGACY?\n")
cat("   Pearson r between legacy profile and averaged-new profile, per timepoint.\n")
cat("   Low r in CM but high r in SH => the NEW CM batch specifically drifted.\n")
cat("=====================================================================\n")
for (ln_ in c("CM","SH")) for (tt in c("control","24h","48h")) {
  a <- paste(ln_,tt,"legacy",sep="_"); b <- paste(ln_,tt,"new",sep="_")
  if (a %in% colnames(avg) & b %in% colnames(avg))
    cat(sprintf("   %-3s %-8s  legacy-vs-new  r = %.3f\n", ln_, tt, cor(avg[,a], avg[,b])))
}

cat("\n=====================================================================\n")
cat("B) WITHIN-BATCH REPRODUCIBILITY: do the two NEW replicates agree?\n")
cat("   Pearson r between new rep1 and new rep2, per condition.\n")
cat("=====================================================================\n")
for (ln_ in c("CM","SH")) for (tt in c("control","24h","48h")) {
  cols <- which(dds$line==ln_ & dds$time==tt & dds$batch=="new")
  if (length(cols)==2) cat(sprintf("   %-3s %-8s  rep1-vs-rep2   r = %.3f\n", ln_, tt, cor(ln[,cols[1]], ln[,cols[2]])))
}

cat("\n=====================================================================\n")
cat("C) SIGNAL-TO-NOISE per library (is the library peak-enriched or background?)\n")
cat("   top1pct = share of signal in the strongest 1%% of peaks (higher = crisper)\n")
cat("   pct_zero = share of consensus peaks with 0 reads (higher = flatter/background)\n")
cat("   near_bg  = share of peaks below 4 normalized reads (weak/near-background)\n")
cat("=====================================================================\n")
cnt <- counts(dds); nrmc <- counts(dds, normalized=TRUE)
tab <- data.frame(
  sample = dds$sample, line = dds$line, time = dds$time, batch = dds$batch,
  top1pct = round(apply(cnt,2,function(x){s<-sort(x,decreasing=TRUE);sum(s[1:ceiling(length(s)*.01)])/sum(s)}),3),
  pct_zero= round(colMeans(cnt==0),3),
  near_bg = round(colMeans(nrmc < 4),3)
) %>% arrange(line, batch, time)
print(tab, row.names = FALSE)

cat("\n=====================================================================\n")
cat("D) BIOLOGICAL VARIANCE vs NOISE inside each new-batch condition\n")
cat("   Median across peaks of |rep1-rep2| on the log2 scale. This is literally\n")
cat("   the replicate disagreement DESeq2 turns into 'dispersion'. Bigger = noisier.\n")
cat("=====================================================================\n")
for (ln_ in c("CM","SH")) for (tt in c("control","24h","48h")) {
  cols <- which(dds$line==ln_ & dds$time==tt & dds$batch=="new")
  if (length(cols)==2) cat(sprintf("   %-3s %-8s  median|log2 rep1 - rep2| = %.3f\n", ln_, tt, median(abs(ln[,cols[1]]-ln[,cols[2]]))))
}
