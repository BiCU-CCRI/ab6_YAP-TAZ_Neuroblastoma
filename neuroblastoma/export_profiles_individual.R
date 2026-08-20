#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Export one profile PDF per target protein, named by the factor, so each file
# is self-identifying (instead of one combined PDF with unlabeled pages).
# Output: results/CnR/profile_<FACTOR>.pdf
#
# Run:  Rscript ~/neuroblastoma/export_profiles_individual.R
# ---------------------------------------------------------------------------
.libPaths(c("/home/abykov/neuroblastoma/R/4.5", .libPaths()))
setwd("/home/abykov/neuroblastoma")

suppressMessages({
  library(DiffBind)
  library(profileplyr)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

res_dir <- path.expand("~/neuroblastoma/results/CnR")

folder <- "data/CnR/consensus_peaks"
prots  <- sub("^dba_(.*)\\.RData$", "\\1",
              list.files(folder, pattern = "^dba_.*\\.RData$"))

.hg38_si <- GenomeInfoDb::seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
sanitize_sites <- function(gr, flank = 5000) {
  if (length(gr) == 0) return(gr)
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  common <- intersect(GenomeInfoDb::seqlevels(gr), GenomeInfoDb::seqlevels(.hg38_si))
  gr <- GenomeInfoDb::keepSeqlevels(gr, common, pruning.mode = "coarse")
  GenomeInfoDb::seqlengths(gr) <- GenomeInfoDb::seqlengths(.hg38_si)[common]
  len <- GenomeInfoDb::seqlengths(gr)[as.character(GenomeInfoDb::seqnames(gr))]
  gr[start(gr) - flank >= 1 & end(gr) + flank <= len]
}

for (nm in prots) {
  cat("\n### Object:", nm, "\n")
  dfbobj <- dba.load(file.path(folder, paste0("dba_", nm)),
                     dir = ".", pre = "", ext = "RData")

  rep <- dba.report(dfbobj, contrast = 2)
  repUP  <- rep[rep$Fold > 2, ]
  repDWN <- rep[rep$Fold < -1, ]
  repUP  <- repUP[order(repUP$Fold, decreasing = TRUE), ]
  repDWN <- repDWN[order(repDWN$Fold, decreasing = FALSE), ]

  repUP  <- sanitize_sites(repUP);  seqlevelsStyle(repUP)  <- "Ensembl"
  repDWN <- sanitize_sites(repDWN); seqlevelsStyle(repDWN) <- "Ensembl"
  rep    <- sanitize_sites(rep);    seqlevelsStyle(rep)    <- "Ensembl"
  rep    <- rep[abs(rep$Fold) > 2, ]

  if (length(repDWN) != 0) {
    profiles <- dba.plotProfile(dfbobj, merge = c(DBA_TISSUE, DBA_REPLICATE),
                                contrast = 2, sites = GRangesList(UP = repUP, DWN = repDWN))
  } else {
    profiles <- dba.plotProfile(dfbobj, merge = c(DBA_TISSUE, DBA_REPLICATE),
                                contrast = 2, sites = rep)
  }

  out_pdf <- file.path(res_dir, paste0("profile_", nm, ".pdf"))
  pdf(file = out_pdf, width = 12, height = 6)
  # Draw the heatmap; add the factor name as the column title so the page is
  # labeled even outside its filename. Fall back to a plain draw if the title
  # passthrough is not accepted.
  drawn <- tryCatch({ dba.plotProfile(profiles, column_title = nm); TRUE },
                    error = function(e) FALSE)
  if (!drawn) dba.plotProfile(profiles)
  dev.off()
  cat("  wrote", out_pdf, "\n")
}
cat("\nAll individual profile PDFs written.\n")
