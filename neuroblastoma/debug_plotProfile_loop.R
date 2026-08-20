#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Faithful reproducer of the Rmd profile loop (02_CnR_analysis.Rmd), now
# writing ONE PDF PER TARGET PROTEIN named by the factor, so each file is
# self-identifying:  results/CnR/profile_<FACTOR>.pdf
#
# Each object is wrapped in tryCatch so we see exactly which object / which
# call fails and with what message.
#
# Run:  Rscript ~/neuroblastoma/debug_plotProfile_loop.R
# ---------------------------------------------------------------------------
.libPaths(c("/home/abykov/neuroblastoma/R/4.5", .libPaths()))
setwd("/home/abykov/neuroblastoma")

suppressMessages({
  library(DiffBind)
  library(profileplyr)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

cat("BiocParallel default:\n"); print(BiocParallel::bpparam())

res_dir <- path.expand("~/neuroblastoma/results/CnR")

# --- load all objects like the Rmd ------------------------------------------
folder <- "data/CnR/consensus_peaks"
prots  <- sub("^dba_(.*)\\.RData$", "\\1",
              list.files(folder, pattern = "^dba_.*\\.RData$"))
DBobj_list <- list()
for (p in prots) {
  DBobj_list[[p]] <- dba.load(file.path(folder, paste0("dba_", p)),
                              dir = ".", pre = "", ext = "RData")
}
cat("Loaded objects:", paste(names(DBobj_list), collapse = ", "), "\n")

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

for (nm in names(DBobj_list)) {
  cat("\n=========================================================\n")
  cat("### Object:", nm, "\n")
  dfbobj <- DBobj_list[[nm]]

  out <- tryCatch({
    rep <- dba.report(dfbobj, contrast = 2)

    repUP  <- rep[rep$Fold > 2, ]
    repDWN <- rep[rep$Fold < -1, ]
    repUP  <- repUP[order(repUP$Fold, decreasing = TRUE), ]
    repDWN <- repDWN[order(repDWN$Fold, decreasing = FALSE), ]

    rep    <- sanitize_sites(rep);    seqlevelsStyle(rep)    <- "Ensembl"
    repUP  <- sanitize_sites(repUP);  seqlevelsStyle(repUP)  <- "Ensembl"
    repDWN <- sanitize_sites(repDWN); seqlevelsStyle(repDWN) <- "Ensembl"

    repList <- GRangesList(UP = repUP, DWN = repDWN)
    rep <- rep[abs(rep$Fold) > 2, ]
    cat(sprintf("  sites: rep=%d  UP=%d  DWN=%d\n",
                length(rep), length(repList$UP), length(repList$DWN)))

    # NOTE: the two batches are stored in the Treatment column (batch1/batch2),
    # so DBA_TREATMENT is added to `merge` to collapse both batches into one
    # profile per condition.
    if (length(repList$DWN) != 0) {
      profiles <- dba.plotProfile(dfbobj,
                                  merge = c(DBA_TISSUE, DBA_REPLICATE, DBA_TREATMENT),
                                  contrast = 2, sites = repList)
    } else {
      profiles <- dba.plotProfile(dfbobj,
                                  merge = c(DBA_TISSUE, DBA_REPLICATE, DBA_TREATMENT),
                                  contrast = 2, sites = rep)
    }
    cat("  profiles built OK\n")

    # --- individual, factor-named PDF ---------------------------------------
    out_pdf <- file.path(res_dir, paste0("profile_", nm, ".pdf"))
    pdf(file = out_pdf, width = 6, height = 12)
    drawn <- tryCatch({ dba.plotProfile(profiles, column_title = nm); TRUE },
                      error = function(e) FALSE)
    if (!drawn) dba.plotProfile(profiles)
    dev.off()
    cat("  wrote", out_pdf, "\n")
    "OK"
  },
  error = function(e) {
    cat("  !!! FAILED:", conditionMessage(e), "\n")
    "FAIL"
  })
  cat("### Result for", nm, ":", out, "\n")
}
cat("\nAll done.\n")
