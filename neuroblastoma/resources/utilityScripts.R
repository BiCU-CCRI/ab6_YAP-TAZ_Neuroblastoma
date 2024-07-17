pheatmap.type <- function(Data, annRow, type=colnames(annRow)[1], doTranspose=FALSE, conditions="Auto", ...){
  ## annRow: A data frame with row names the same as row names of Data.
  ## type: The column name of annRow representing two or more conditions.
  ## This function first performs hierarchical clustering on samples
  ## (rows of Data) within each condition.
  ##^Then plots a heatmap without further clustering of rows. 
  ## ... are passed to pheatmap function.
  res <- list()
  annRow <- annRow[, type, drop=FALSE]
  ## QC:
  if(is.null(rownames(annRow)))
    stop("annrow must have row names!")
  if(any(! rownames(annRow) %in% rownames(Data)))
    stop("annRow has rows that are not present in rows of Data!")
  ## Put all samples in the same condition together.
  annRow <- annRow[order(annRow[, 1]), , drop=FALSE]
  samplesOriginalOrder <- rownames(Data)
  ##Data <- Data[rownames(annRow), , drop=FALSE]
  if(conditions[1]=="Auto")
    conditions <- unique(as.character(annRow[, 1]))
  if(any(!conditions %in% unique(as.character(annRow[, 1])))){
    warning("Some of the conditions are not in annRow.")
    conditions <- intersect(conditions, unique(as.character(annRow[, 1])))
  }
  pheatmapS <- list()
  dataPlot <- c()
  ann1 <- c()
  for(cond in conditions){
    condSamples <- rownames(annRow)[which(annRow==cond)]  
    if(length(condSamples)>1){
      pa <- pheatmap(Data[condSamples, , drop=FALSE], cluster_cols=FALSE, silent=TRUE)
      pheatmapS[[as.character(cond)]] <- pa
      o2 <- pa$tree_row$order
    } else {
      o2 <- 1
    }
    dataPlot <- rbind(dataPlot, Data[condSamples[o2], , drop=FALSE])
    ann1 <- rbind(ann1, annRow[condSamples[o2], , drop=FALSE])
  }
  if(!doTranspose){
    pAll <- pheatmap(dataPlot, annotation_row=ann1, cluster_rows=FALSE, ...)
  } else { ## Transpose
    pAll <- pheatmap(t(dataPlot), annotation_col=ann1, cluster_cols=FALSE, ...)
  }
  res[["pheatmapS"]] <- pheatmapS
  res[["pheat"]] <- pAll
  res[["ordering"]] <- match(rownames(dataPlot), samplesOriginalOrder)
  res[["annRowAll"]] <- ann1
  invisible(res)
}



generatePCA <- function(transf_object=NULL, cond_interest_varPart=NULL, color_variable=NULL, shape_variable=NULL, ntop_genes=500){
  
  library(DESeq2)
  library(ggplot2)
  
  #transf_object_counts <- assay(transf_object)
  #ntop_genes=nrow(transf_object_counts) 
  pcaData <- DESeq2::plotPCA(transf_object, intgroup=cond_interest_varPart, returnData=TRUE, ntop=ntop_genes) 
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ggplot(pcaData, aes(PC1, PC2, color=!!sym(color_variable), shape=!!sym(shape_variable))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #+ coord_fixed()
    theme_bw() 
  
}


meanExprsPerGroup <- function(dds_object=NULL, 
                              #condition_test=NULL,
                              variable=NULL,
                              cond_numerator=NULL, 
                              cond_denominator=NULL){
  require(DESeq2)
  require(dplyr)
  # extracts expression matrix (all samples) for numerator and denominator
  # calculates mean expression per group (numerator, denominator)
  # [ ] add check that there is filenames and condition column!
  # Normalized counts and means per group
  # creating only subset for actual comparison
  # may need to change condition to some other variable
  # res_extract <- condition_test
  # cond_numerator <- gsub(pattern = paste0("(", variable,"_)(.+)(_vs_.+)"), replacement = "\\2", res_extract)
  # cond_denominator <- gsub(pattern = paste0("(", variable,"_)(.+_vs_)(.+)"), replacement = "\\3", res_extract)
  
  new_sample_names <- as.data.frame(colData(dds_object)) 
  
  if("filenames" %in% names(new_sample_names)){
    new_sample_names <- new_sample_names %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
    
  } else {
    new_sample_names <- new_sample_names %>%
      tibble::rownames_to_column(., var = "filenames") %>%
      dplyr::select(filenames, tidyselect::all_of(variable)) 
  }
  
  # check if variable is factor otherwise arranging in the next step will be "random"
  #if (is.factor()) {
  #  
  #}
  
  # extract filenames for each of the conditions and pivot_table
  new_sample_names <- new_sample_names %>%
    dplyr::select(filenames, !!as.name(variable)) %>% # selecting filenames and variable of interest (e.g.)
    dplyr::filter(!!as.name(variable) == cond_numerator | !!as.name(variable) == cond_denominator) %>% # filtering to keep only numerator and denominator samples
    dplyr::transmute(filenames,
                     denom_num_extract = factor(!!as.name(variable), levels=c(cond_denominator, cond_numerator))) %>% # re-factoring denominator, then numerator (but this should have been done in cond_data already)
    dplyr::arrange(., denom_num_extract) %>% # convert the strings to names with as.name and !! unquote (bang-bang) !!as.name(variable)
    # reorder variable according to numerator, denominator; so the expression output is in correct order - this assumes previous correct ordering
    dplyr::mutate(new_name = paste(denom_num_extract, filenames, sep="-")) # creating new column with variable of interest and filename
  
  
  normalized_counts <- NULL # just to make sure it does not exist from previous run
  normalized_counts <- data.frame(counts(dds_object, normalized = TRUE))
  # extracting subset
  normalized_counts <- normalized_counts %>%
    dplyr::select(new_sample_names$filenames) # keeping only conditions that are being compared and following previous order
  
  colnames(normalized_counts) <- new_sample_names$new_name
  normalized_counts <- normalized_counts %>%
    tibble::rownames_to_column("peak_id")
  
  normalized_counts_AddedMean <- normalized_counts %>%
    dplyr::mutate(., 
                  MeanExpr_denominator = rowMeans(dplyr::select(., matches(paste0(cond_denominator,"-"))), na.rm = TRUE),
                  MeanExpr_numerator = rowMeans(dplyr::select(., matches(paste0(cond_numerator,"-"))), na.rm = TRUE)) # more robust regex?!
  
  # rename mean_numerator, mean_denominator in the final column
  # use rename_all! rename(new_sample_names$new_sample_names)
  # rename colnames to condition_Sample number name
  # calculate mean across normalized counts
  
  return(normalized_counts_AddedMean)
}


extract_results_DDS <- function(dds_object = NULL,
                                coeff_name = NULL,
                                cond_numerator = NULL,
                                cond_denominator = NULL,
                                cond_variable = NULL,
                                padj_cutoff = 0.5,
                                log2FC_cutoff = 0.58){
  
  #This is the start of the PR function that I re-purposed.
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          name=coeff_name,   
                                          parallel = TRUE, 
                                          alpha = padj_cutoff)
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,   
                                 res=res_table_unshrunken, 
                                 type = "apeglm")
  row.names(res_table) <-   dds_object@rowRanges@elementMetadata@listData[["peak_id"]] #DATASET   SPECIFIC!
  
  normalized_counts_AddedMean <-   meanExprsPerGroup(dds_object=dds_object,
                                                     cond_numerator=cond_numerator,
                                                     cond_denominator=cond_denominator,
                                                     variable = cond_variable)
  normalized_counts_AddedMean$peak_id <-   dds_object@rowRanges@elementMetadata@listData[["peak_id"]]
  
  ensemblAnnot <- as.data.frame(dds_object@rowRanges@elementMetadata@listData)
  ######
  #Operations with table, merging with metadata.
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("peak_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, 
                                    -2^(abs(log2FoldChange)), 
                                    2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "peak_id") %>% # adding annotationl   entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "peak_id") %>% # adding   normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(peak_id, gencode_chr, gencode_start, gencode_end, 
                  sEnh_type, YAP_TAZ, is_promoter_3kb, gene_category_our_RNAseq, gene_category_GROEN, 
                  gencode_gene_name, gencode_characterization, homer_Annotation, 
                  homer_Distance.to.TSS, homer_Nearest.PromoterID, homer_Nearest.Ensembl, 
                  baseMean.y, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, pvalue, padj,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns   not mentioned; arts_with(paste0(cond_denominator,"_S") may not work if other   references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs =   funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find   better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs =   funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design =   paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes =   nrow(results_data_annot_signif),
                                        signif_genes_UP =   sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN =   sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) &   (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  dds_consens_results_list<-list(results_signif=results_data_annot_signif  , de_details=temp_results_summary_df, results_all=results_data_annot)
  return(dds_consens_results_list)
  #END OF PR FUNCTION
  
}




plotVolcano <- function(dds_results_obj=NULL, 
                        genes_of_interest=NULL, 
                        plot_title=NULL, 
                        log2FC_cutoff=0.58, 
                        padj_cutoff=0.05)
{ 
  #browser()
  results_data_annot_forPlot <- dds_results_obj %>% 
    dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "DOWN"
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE,
                                                 levels = c("DOWN", "NO", "UP"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point(data = subset(results_data_annot_forPlot, signif_DE == 'NO'), 
               aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point(data = subset(results_data_annot_forPlot, signif_DE == 'DOWN'), 
               aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    geom_point(data = subset(results_data_annot_forPlot, signif_DE == 'UP'), 
               aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
    ggrepel::geom_label_repel(data = . %>% 
                                filter(gencode_gene_name %in% genes_of_interest) %>%
                                filter(abs(log2FoldChange) > log2FC_cutoff) %>%
                                filter(padj < padj_cutoff), 
                              aes(label = gencode_gene_name),
                              show.legend = FALSE,
                              box.padding = 0.5,
                              segment.color ="black",
                              max.overlaps = Inf
                              ) +
    scale_color_manual(values=c(NO = "grey", DOWN="navy", UP="firebrick3")) +
    theme_bw(base_size = 14) +
    labs(x = "log2FC") + 
    ggtitle(plot_title)
  
  return(signif_volcanoPlot)
}







extract_results_DDS_HIC <- function(dds_object = NULL,
                                    coeff_name = NULL,
                                    cond_numerator = NULL,
                                    cond_denominator = NULL,
                                    cond_variable = NULL,
                                    padj_cutoff = 0.5,
                                    log2FC_cutoff = 0.58){
  
  #This is the start of the PR function that I re-purposed.
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          name=coeff_name,   
                                          parallel = TRUE, 
                                          alpha = padj_cutoff)
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,   
                                 res=res_table_unshrunken, 
                                 type = "apeglm")
  row.names(res_table) <-   dds_object@rowRanges@elementMetadata@listData[["peak_id"]] #DATASET   SPECIFIC!
  
  normalized_counts_AddedMean <-   meanExprsPerGroup(dds_object=dds_object,
                                                     cond_numerator=cond_numerator,
                                                     cond_denominator=cond_denominator,
                                                     variable = cond_variable)
  normalized_counts_AddedMean$peak_id <-   dds_object@rowRanges@elementMetadata@listData[["peak_id"]]
  
  ensemblAnnot <- as.data.frame(dds_object@rowRanges@elementMetadata@listData)
  ######
  #Operations with table, merging with metadata.
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("peak_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, 
                                    -2^(abs(log2FoldChange)), 
                                    2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "peak_id") %>% # adding annotationl   entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "peak_id") %>% # adding   normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(peak_id, gencode_chr, gencode_start, gencode_end, 
                  sEnh_type, YAP_TAZ, is_promoter_3kb,  gene_category_our_RNAseq, gene_category_GROEN, 
                  gencode_gene_name, gencode_characterization, homer_Annotation, ENH, 
                  homer_Distance.to.TSS, homer_Nearest.PromoterID, homer_Nearest.Ensembl, 
                  baseMean.y, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, pvalue, padj,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns   not mentioned; arts_with(paste0(cond_denominator,"_S") may not work if other   references names dont follow with _S
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs =   funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find   better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs =   funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
                                        design =   paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes =   nrow(results_data_annot_signif),
                                        signif_genes_UP =   sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN =   sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) &   (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  
  dds_consens_results_list<-list(results_signif=results_data_annot_signif  , de_details=temp_results_summary_df, results_all=results_data_annot)
  return(dds_consens_results_list)
  #END OF PR FUNCTION
  
}




import:::here(gseaScores, .from = DOSE)
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin=0
  df$ymax=0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

tableGrob2 <- function(d, p = NULL) {
  d <- d[order(rownames(d)),]
  tp <- gridExtra::tableGrob(d) # from gridExtra
  if (is.null(p)) {
    return(tp)
  }
  pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] = grid::gpar(col = pcol[i])
  }
  return(tp)
}

gseaplot3 <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                       rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                       ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "NES", "p.adjust")]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    for (i in seq_len(ncol(pd))) {
      pd[, i] <- format(pd[, i], digits = 4)
    }
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.9))
  }
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_text())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  aplot::gglist(gglist = plotlist, ncol = 1, heights = rel_heights)
}








chipEnrichAndExport <- function(peaks, peaksName, locusdef, res_dir, TF_name, genesets, genesets_name ) {
  
  results <- chipenrich(
    peaks = peaks,
    genome = "hg38",
    genesets = genesets,
    locusdef = locusdef,
    qc_plots = TRUE,
    out_name = NULL,
    n_cores = 1,
    # randomization = 'complete',
    max_geneset_size = 5000
  )
  
  XLSX_OUT <- createWorkbook()
  addWorksheet(XLSX_OUT, "results")
  addWorksheet(XLSX_OUT, "peaks")
  addWorksheet(XLSX_OUT, "peaks_per_gene")
  addWorksheet(XLSX_OUT, "opts")
  
  writeData(XLSX_OUT, x = results$results, sheet = "results")
  writeData(XLSX_OUT, x = results$peaks, sheet = "peaks")
  writeData(XLSX_OUT, x = results$peaks_per_gene, sheet = "peaks_per_gene")
  writeData(XLSX_OUT, x = results$opts, sheet = "opts")
  
  saveWorkbook(XLSX_OUT,
               file.path(res_dir, paste0("Enricher_rslts_",  TF_name, "_", genesets_name, "_", locusdef, "_", peaksName, "_peaks.xlsx")),
               overwrite = TRUE)

}

