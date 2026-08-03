#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Build CUT&RUN consensus peaks with DiffBind from MACS2 narrowPeak files.
#
# This mirrors the LEGACY recipe used to create the original dba_<Factor>.RData
# objects, so the regenerated objects are structurally identical:
#   dba() -> dba.count() -> dba.normalize() -> dba.contrast() -> dba.analyze()
# giving each object the design / norm / contrasts / DESeq2 / blacklist /
# greylist slots the notebook (02_CnR_analysis.Rmd) depends on.
#
# Per factor (H3K27ac, H3K4me1, Jun, TAZ, YAP) it:
#   1. assembles a DiffBind sample sheet (ChIP bam + narrowPeak + summits + IgG),
#   2. counts reads over the default consensus,
#   3. normalizes,
#   4. builds contrasts on design ~Treatment + Tissue + Condition, where
#      Treatment = batch (1 vs 2); contrast 1 = Tissue, contrast 2 = Condition
#      (MES vs ADRN) -- the `contrast = 2` the notebook uses,
#   5. runs the differential analysis (with blacklist + greylist, as legacy did),
#   6. saves `dba_<Factor>.RData` and a `<Factor>_consenus_peaks.bed` for HOMER.
#
# NOTE: the legacy pipeline additionally annotated the consensus BED with HOMER
# (annotatePeaks.pl, run OUTSIDE R) and merged the result back into
# dfbobj$DESeq2$DEdata mcols. That external step is not reproduced here; the
# objects written are equivalent to the legacy "_not_annotated" stage but saved
# under the final `dba_<Factor>` name the notebook loads.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DiffBind)
  library(stringr)
  library(GreyListChIP)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

# --- greylist fix (from 02_CnR_analysis.Rmd chunk 1) -----------------------
# The default greylist step chokes when the BAM header contains scaffold
# chromosomes not in the karyotype; restrict to the intersection.
pv.countGreylistEdited <- function(bamfile, pv, ktype) {
  gl <- new("GreyList", karyotype = ktype[intersect(pv$chrmap, names(ktype)), ])
  gl <- GreyListChIP::countReads(gl, bamfile)
  return(gl)
}
environment(pv.countGreylistEdited) <- asNamespace("DiffBind")
assignInNamespace("pv.countGreylist", pv.countGreylistEdited, ns = "DiffBind")

# ---- paths & options (edit if you move things) ----------------------------
peaks_dir  <- "/home/abykov/neuroblastoma/data/CnR/consensus_peaks"
bam_dir    <- "/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/raw_data/CnR/outdir/02_alignment/bowtie2/target/markdup"
out_dir    <- peaks_dir                    # where dba_<Factor> + BEDs are written

design_formula <- "~Treatment + Tissue + Condition"  # Treatment carries batch (1 vs 2)
apply_blacklist_greylist <- TRUE           # legacy default; set FALSE for a quick run without masking

bam_suffix <- ".target.markdup.sorted.bam"
marks      <- c("H3K27ac", "H3K4me1", "Jun", "TAZ", "YAP")

# ---- 1. assemble the sample sheet from the narrowPeak files ---------------
peak_files   <- list.files(peaks_dir, pattern = "\\.macs2_peaks\\.narrowPeak$",
                          full.names = TRUE)
summit_files <- sub("\\.macs2_peaks\\.narrowPeak$", ".macs2_summits.bed", peak_files)
base         <- sub("\\.macs2_peaks\\.narrowPeak$", "", basename(peak_files))  # e.g. SK-N-SH-M_H3K27ac_R1

rep_num  <- as.integer(str_match(base, "_R(\\d+)$")[, 2])
peakbase <- sub("_R\\d+$", "", base)                                          # e.g. SK-N-SH-M_H3K27ac

Factor    <- str_extract(peakbase, paste(marks, collapse = "|"))
Tissue    <- str_extract(peakbase, "SK-N-SH|CLB-M[Aa]")
Tissue    <- ifelse(grepl("^CLB", Tissue), "CLB-Ma", Tissue)                  # normalise CLB-MA -> CLB-Ma
state     <- str_match(peakbase, "^(?:SK-N-SH|CLB-M[Aa])-(A|M|ADR|MES)[-_]")[, 2]
Condition <- ifelse(state %in% c("A", "ADR"), "ADRN", "MES")                  # A/ADR -> ADRN ; M/MES -> MES

# Batch: original A/M naming = batch 1 ; ADR/MES naming = batch 2.
# Stored in DiffBind's Treatment metadata slot so it can enter the design formula
# (referenced as `Treatment`). Batch is crossed with Tissue and Condition here,
# so it is estimable as a covariate.
Batch     <- ifelse(state %in% c("A", "M"), "batch1", "batch2")

# control (IgG) of the same cell-line + state: swap the mark token for IgG
ctrl_group <- sub(paste0("(", paste(marks, collapse = "|"), ")$"), "IgG", peakbase)

# replicate-matched control bam, falling back to _R1 when a replicate is missing
# (legacy did the same for CLB-Ma-M_IgG_R2 -> _R1)
ctrl_bam <- file.path(bam_dir, paste0(ctrl_group, "_R", rep_num, bam_suffix))
missing  <- !file.exists(ctrl_bam)
ctrl_bam[missing] <- file.path(bam_dir, paste0(ctrl_group[missing], "_R1", bam_suffix))

samples <- data.frame(
  SampleID    = base,
  Tissue      = Tissue,
  Factor      = Factor,
  Condition   = Condition,
  Treatment   = Batch,          # batch (1 vs 2) as a modeled covariate
  Replicate   = rep_num,
  bamReads    = file.path(bam_dir, paste0(base, bam_suffix)),
  ControlID   = paste0(ctrl_group, "_R", rep_num),
  bamControl  = ctrl_bam,
  Peaks       = peak_files,
  PeakCaller  = "narrow",
  PeakSummits = summit_files,
  stringsAsFactors = FALSE
)
samples$Peaks_count <- sapply(peak_files, function(f) length(readLines(f, warn = FALSE)))

# sanity: every referenced file must exist
stopifnot(all(file.exists(samples$bamReads)))
stopifnot(all(file.exists(samples$bamControl)))
stopifnot(all(file.exists(samples$Peaks)))
stopifnot(all(file.exists(samples$PeakSummits)))
message(sprintf("Assembled sample sheet: %d samples across %d factors",
                nrow(samples), length(unique(samples$Factor))))
write.csv(samples, file.path(out_dir, "CnR_diffbind_samplesheet.csv"), row.names = FALSE)

# ---- 2-6. per-factor: count -> normalize -> contrast -> analyze -> save ----
DBobj_list <- list()
for (mark in unique(samples$Factor)) {
  message("=== ", mark, " ===")
  s <- samples[samples$Factor == mark, ]

  dbo <- dba(sampleSheet = s)          # loads peaksets
  dbo <- dba.count(dbo)                # default consensus/minOverlap, like legacy
  dbo <- dba.normalize(dbo)            # populates $norm

  # Explicit design mode so batch (Treatment) is controlled for while keeping the
  # contrast order the notebook expects. Each dba.contrast(contrast=) call appends
  # one contrast -> contrast 1 = Tissue, contrast 2 = Condition (MES vs ADRN).
  # The explicit c("Condition","MES","ADRN") sets direction: Fold > 0 = higher in
  # MES (matches the notebook's select_up_or_down_regulated = "up" == MES).
  dbo <- dba.contrast(dbo, design = design_formula,
                     contrast = c("Tissue", "SK-N-SH", "CLB-Ma"))    # contrast 1
  dbo <- dba.contrast(dbo, contrast = c("Condition", "MES", "ADRN")) # contrast 2

  if (apply_blacklist_greylist) {
    dbo <- dba.analyze(dbo)                                     # default: blacklist + greylist
  } else {
    dbo <- dba.analyze(dbo, bBlacklist = FALSE, bGreylist = FALSE)
  }

  message("  contrasts (verify contrast 2 = Condition MES vs ADRN):")
  print(dba.show(dbo, bContrasts = TRUE))
  DBobj_list[[mark]] <- dbo

  # save the DiffBind object under the name 02_CnR_analysis.Rmd expects
  dba.save(dbo, file = paste0("dba_", mark), dir = out_dir, pre = "", ext = "RData")

  # export consensus peaks in the legacy BED layout for HOMER annotatePeaks.pl:
  # columns = peakID, chr, start, end, "."  (no header)
  cons <- as.data.frame(dba.peakset(dbo, bRetrieve = TRUE))[, 1:3]
  cons <- cbind(rownames(cons), cons, ".")
  cons[[2]] <- ifelse(grepl("^chr", cons[[2]]), cons[[2]], paste0("chr", cons[[2]]))
  write.table(cons, file.path(out_dir, paste0(mark, "_consenus_peaks.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  message(sprintf("  consensus peaks: %d", nrow(cons)))
}

saveRDS(DBobj_list, file = "~/neuroblastoma/results/CnR/DBobj_list.RDS")

message("Done. Wrote dba_<Factor>.RData and <Factor>_consenus_peaks.bed to ", out_dir)
