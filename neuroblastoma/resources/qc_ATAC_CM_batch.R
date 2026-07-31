#!/usr/bin/env Rscript
# QC diagnostics for the ATAC-seq k975 drug-treatment analysis.
# Goal: understand why the new CM-line samples behave oddly and why nothing
# passes DESeq2 significance filters in 03_ATAC-seq_k975_inhibition.Rmd.
#
# Run:  R_LIBS_USER=/home/abykov/R/4.5 Rscript resources/qc_ATAC_CM_batch.R
.libPaths(c("/home/abykov/R/4.5", .libPaths()))

suppressMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(dplyr)
  library(stringr)
})

out_dir <- "~/neuroblastoma/results/ATAC-seq/ATAC-seq_drug_treatment/QC"
out_dir <- path.expand(out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_con <- file(file.path(out_dir, "qc_report.txt"), open = "wt")
say <- function(...) { cat(..., "\n"); cat(..., "\n", file = log_con) }

# ---------------------------------------------------------------------------
# 1. Load and subset exactly as the Rmd does
# ---------------------------------------------------------------------------
se <- readRDS("~/neuroblastoma/data/ATACseq/consensus_peaks.mLb.clN.rds")

samples_to_subset <- c(
  "CM_control_ATAC_S165169_REP1", "CM_24h_ATAC_S165170_REP1", "CM_48h_ATAC_S165167_REP1",
  "SH_control_ATAC_S165166_REP1", "SH_24h_ATAC_S165164_REP1", "SH_48h_ATAC_S165165_REP1",
  "SH_DMSO_ctrl_rep1_REP1", "SH_DMSO_ctrl_rep2_REP1",
  "SH_24h_rep1_REP1", "SH_24h_rep2_REP1", "SH_48h_rep1_REP1", "SH_48h_rep2_REP1",
  "CM_DMSO_ctrl_rep1_REP1", "CM_DMSO_ctrl_rep2_REP1",
  "CM_24h_rep1_REP1", "CM_24h_rep2_REP1", "CM_48h_rep1_REP1", "CM_48h_rep2_REP1"
)

dds <- se[, colData(se)$sample %in% samples_to_subset]

colData(dds)$cell_line <- factor(case_when(
  str_detect(dds$sample, "^CM") ~ "CLBMa",
  str_detect(dds$sample, "^SH") ~ "SK_N_SH",
  TRUE ~ as.character(dds$sample)
))
colData(dds)$time <- factor(case_when(
  str_detect(dds$sample, "24h") ~ "24h",
  str_detect(dds$sample, "48h") ~ "48h",
  str_detect(dds$sample, "DMSO") ~ "control",
  str_detect(dds$sample, "control") ~ "control",
  TRUE ~ NA_character_
))
colData(dds)$time <- relevel(colData(dds)$time, ref = "control")
colData(dds)$batch <- factor(if_else(str_detect(dds$sample, "_ATAC_S\\d+"), "legacy", "new"))
colData(dds)$grp <- paste(dds$cell_line, dds$time, dds$batch, sep = "_")

meta <- as.data.frame(colData(dds))
say("=========================================================")
say("SAMPLE / DESIGN TABLE")
say("=========================================================")
print(table(meta$cell_line, meta$time, meta$batch))
capture.output(print(table(meta$cell_line, meta$time, meta$batch)), file = log_con)

# ---------------------------------------------------------------------------
# 2. Library size / total signal per sample  (peak-count depth)
# ---------------------------------------------------------------------------
cnt <- assay(dds, "counts")
libsize <- colSums(cnt)
meta$libsize <- libsize[rownames(meta)]
# fraction of the total signal captured by the top 1% strongest peaks per sample
top_frac <- apply(cnt, 2, function(x) {
  s <- sort(x, decreasing = TRUE)
  sum(s[seq_len(ceiling(length(s) * 0.01))]) / sum(s)
})
meta$top1pct_frac <- top_frac[rownames(meta)]
# fraction of peaks with zero counts
meta$frac_zero_peaks <- colMeans(cnt == 0)[rownames(meta)]

say("\n=========================================================")
say("PER-SAMPLE DEPTH & SIGNAL CONCENTRATION")
say("  libsize        = total reads in consensus peaks")
say("  top1pct_frac   = fraction of signal in top 1% peaks (high = concentrated/good)")
say("  frac_zero_peaks= fraction of peaks with 0 counts (high = sparse/noisy)")
say("=========================================================")
qc_tab <- meta %>%
  dplyr::select(sample, cell_line, time, batch, libsize, top1pct_frac, frac_zero_peaks) %>%
  arrange(cell_line, batch, time)
print(qc_tab, row.names = FALSE)
capture.output(print(qc_tab, row.names = FALSE), file = log_con)
write.csv(qc_tab, file.path(out_dir, "per_sample_depth_qc.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Filter + VST, then PCA and sample correlation
# ---------------------------------------------------------------------------
design(dds) <- ~ batch + cell_line + time
keep <- rowSums(fpm(dds, robust = TRUE) > 1) >= 2
say("\nPeaks kept after filter (fpm>1 in >=2 samples):", sum(keep), "of", length(keep))
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)

say("\nDESeq2 size factors (median-of-ratios). Extreme values => composition problem:")
sf <- sizeFactors(dds)
sf_tab <- data.frame(sample = names(sf), sizeFactor = round(sf, 3),
                     cell_line = meta[names(sf), "cell_line"],
                     batch = meta[names(sf), "batch"])[order(sf), ]
print(sf_tab, row.names = FALSE)
capture.output(print(sf_tab, row.names = FALSE), file = log_con)

vsd <- vst(dds, blind = TRUE)
vm <- assay(vsd)

# ---- Sample-sample Pearson correlation on most variable peaks ----
rv <- matrixStats::rowVars(vm)
top <- head(order(rv, decreasing = TRUE), 2000)
cormat <- cor(vm[top, ])

# For each sample, its nearest neighbour by correlation (excluding itself)
nn <- apply(cormat, 1, function(x) {
  x2 <- x; x2[which.max(x2)] <- -Inf   # drop self (==1)
  names(which.max(x2))
})
meta$nearest_neighbour <- nn[rownames(meta)]
meta$nn_expected_same_line <- str_sub(meta$sample, 1, 2) == str_sub(meta$nearest_neighbour, 1, 2)

say("\n=========================================================")
say("NEAREST NEIGHBOUR BY CORRELATION (top 2000 variable peaks)")
say("If a CM sample's nearest neighbour is an SH sample -> swap/contamination flag")
say("=========================================================")
nn_tab <- meta %>%
  dplyr::select(sample, cell_line, time, batch, nearest_neighbour, nn_expected_same_line) %>%
  arrange(cell_line, batch, time)
print(nn_tab, row.names = FALSE)
capture.output(print(nn_tab, row.names = FALSE), file = log_con)

# ---- Within-group replicate correlation (new batch has 2 reps) ----
say("\n=========================================================")
say("REPLICATE CORRELATION within each cell_line x time x batch group")
say("Low values in CM-new groups => noisy/irreproducible libraries")
say("=========================================================")
for (g in unique(meta$grp)) {
  s <- rownames(meta)[meta$grp == g]
  if (length(s) >= 2) {
    cc <- cormat[s, s]
    say(sprintf("  %-28s reps=%d  mean_r=%.3f", g, length(s),
                mean(cc[lower.tri(cc)])))
  }
}

# ---------------------------------------------------------------------------
# 4. Plots: PCA (all + CM only) and correlation heatmap
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(pheatmap) })

pcaData <- plotPCA(vsd, intgroup = c("cell_line", "time", "batch"), returnData = TRUE, ntop = 2000)
pv <- round(100 * attr(pcaData, "percentVar"))
p_all <- ggplot(pcaData, aes(PC1, PC2, color = time, shape = batch)) +
  geom_point(size = 3) +
  facet_wrap(~cell_line) +
  labs(x = paste0("PC1: ", pv[1], "%"), y = paste0("PC2: ", pv[2], "%"),
       title = "ATAC drug-treatment PCA (blind VST, top2000)") +
  theme_bw()
ggsave(file.path(out_dir, "PCA_by_cellline_time_batch.png"), p_all, width = 24, height = 12, units = "cm")

ann <- meta[, c("cell_line", "time", "batch")]
pheatmap(cormat,
         annotation_col = ann, annotation_row = ann,
         main = "Sample-sample Pearson correlation (top2000 variable peaks)",
         fontsize = 6,
         filename = file.path(out_dir, "sample_correlation_heatmap.png"),
         width = 11, height = 10)

# ---- CM-only re-analysis: does the drug signal exist within CM alone? ----
say("\n=========================================================")
say("CM-ONLY DESeq2 (design ~ batch + time) -- does ANY signal survive")
say("when SH is removed and CM is modelled on its own?")
say("=========================================================")
run_line <- function(line_lab) {
  sub <- dds[, dds$cell_line == line_lab]
  colData(sub) <- droplevels(colData(sub))
  design(sub) <- ~ batch + time
  sub <- tryCatch(DESeq(sub, quiet = TRUE), error = function(e) { say("  DESeq error:", conditionMessage(e)); NULL })
  if (is.null(sub)) return(invisible())
  for (cn in c("48h", "24h")) {
    r <- results(sub, contrast = c("time", cn, "control"))
    say(sprintf("  %s  %s vs control:  padj<0.05 = %d ; padj<0.05 & |LFC|>1 = %d ; min padj = %.2e",
                line_lab, cn, sum(r$padj < 0.05, na.rm = TRUE),
                sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE),
                min(r$padj, na.rm = TRUE)))
  }
  say(sprintf("  %s dispersion (median gene-est): %.3f", line_lab,
              median(mcols(sub)$dispGeneEst, na.rm = TRUE)))
}
run_line("CLBMa")
run_line("SK_N_SH")

# ---- pooled model dispersion for reference ----
say("\nPooled model (~batch+cell_line+time) dispersion for comparison:")
dds_pool <- DESeq(dds, quiet = TRUE)
say(sprintf("  median gene-est dispersion (pooled): %.3f",
            median(mcols(dds_pool)$dispGeneEst, na.rm = TRUE)))
for (cn in c("48h", "24h")) {
  r <- results(dds_pool, contrast = c("time", cn, "control"))
  say(sprintf("  pooled %s vs control: padj<0.05 = %d ; &|LFC|>1 = %d ; min padj = %.2e",
              cn, sum(r$padj < 0.05, na.rm = TRUE),
              sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm = TRUE),
              min(r$padj, na.rm = TRUE)))
}

say("\nDone. Outputs in:", out_dir)
close(log_con)
