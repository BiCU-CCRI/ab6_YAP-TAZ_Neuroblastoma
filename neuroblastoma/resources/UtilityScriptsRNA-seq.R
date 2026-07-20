#' generateEnsemblAnnotation: Generates ensembl annotation file for genes in RNA-seq experiment
#' 
#' \code{generateEnsemblAnnotation} constructs, for example, DESeq2 dds object from Salmon count files, condition table and design
#' @param ensembl_ids ensembl_ids of genes in an experiment
#' @param host ensembl website address (e.g. http://www.ensembl.org)
#' @param dataset e.g. mmusculus_gene_ensembl
#' @param version version of dataset

generateEnsemblAnnotation <- function(ensembl_ids=NULL, host="http://www.ensembl.org", dataset="hsapiens_gene_ensembl", version="Ensembl Genes 109"){
  require("biomaRt") # Gene annotation using ensembl database; carefully choose same version that was used for alignment
  require("RCurl") # proxy settings for biomaRt  
  #options(RCurlOptions = list(proxy="specify-proxy-address",http.version=HTTP_VERSION_1_0)) 
  require("dplyr")
  
  ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    dataset=dataset, 
                    host=host, 
                    port = 80, 
                    verbose =T, 
                    version=version)
  
  mart <- useDataset(dataset = dataset, ensembl)
  
  # check the available "filters" - things you can filter for
  #listFilters(ensembl) %>% filter(str_detect(name, "ensembl"))
  filterType <- "ensembl_gene_id"
  # check the available "attributes" - things you can retreive
  #listAttributes(ensembl) %>%  head(20)
  # Set the list of attributes
  
  # temporary fix to annotate also human!
  if(dataset=="hsapiens_gene_ensembl") {
    cat("Annotating: ", dataset, "\n")
    # entrezgene
    attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol', 'description', 'gene_biotype', 
                        'chromosome_name', 'start_position', 'end_position', 'strand')
  } else {
    attributeNames <- c('ensembl_gene_id', 'entrezgene_id', 'mgi_symbol', 'mgi_description', 'gene_biotype',
                        'chromosome_name', 'start_position', 'end_position', 'strand')
  }

  if (!is.null(ensembl_ids) & is.character(ensembl_ids)) {
    print(paste0("Generating annotation file for genes: ", length(ensembl_ids)))
    print(paste0("dataset: ", dataset))
    print(paste0("version: ", version))
    print(paste0("host: ", host))
    
    # Adding gene names position, etc.
    ensembl_ids <- ensembl_ids
    
    # uniprot swissprot
    # add attribute: mgi_symbol
    ensemblAnnot <- getBM(
      filters = filterType, 
      attributes = attributeNames,
      values = ensembl_ids,
      mart = mart)
    
    print("Removing duplicates from ensemblAnnot file...this may take a moment!")
    ensemblAnnot_dupRemove <- ensemblAnnot %>%
      dplyr::rename(ensembl_id = ensembl_gene_id) %>%
      dplyr::group_by(ensembl_id) %>%
      dplyr::summarise_all(., .funs = function(x) {paste(unique(x), collapse = ";")}) # removing duplicates in ensembl_id by grouping; collapsing unique by ";" 
    ensemblAnnot_dupRemove <- as.data.frame(ensemblAnnot_dupRemove)
    rownames(ensemblAnnot_dupRemove) <- ensemblAnnot_dupRemove$ensembl_id
    
    # converting "NA" strings to NA
    ensemblAnnot_dupRemove$entrezgene <- ifelse(ensemblAnnot_dupRemove$entrezgene == "NA", 
                                                NA_character_, 
                                                ensemblAnnot_dupRemove$entrezgene)
    
    return(ensemblAnnot_dupRemove)
    
  } else {
    print("ensembl_gene_id needs to be specified as a character vector.")
  }
}





filterDatasets <- function(dds_object = NULL, 
                           abs_filt = TRUE, 
                           abs_filt_samples = 3, 
                           relat_filt = 0.2) {
  require("DESeq2") # and dependencies to work with dds object
  # filtering samples based on fpm
  # abs_filt = TRUE, filtering on absolute number of samples defined in abs_filt_samples
  # abs_filt = FALSE, filtering on relative number of samples
  
  if (!is.null(dds_object)) {
    cat("Original dds object samples: ", ncol(dds_object), " genes: ", nrow(dds_object), "\n")
    
    if (abs_filt == TRUE) {
      # at least 1 fpm in at least n-number of samples
      # by default considering smallest condition (cluster) with 3 samples
      
      cat("Minimum number of samples with expression:", abs_filt_samples, "\n")
      keep_genes_idx <- rowSums(DESeq2::fpm(dds_object, robust = TRUE) > 1) >= abs_filt_samples
      cat("Number of filtered genes:", sum(keep_genes_idx == FALSE), "\n")
      
      dds_object_filt <- dds_object[keep_genes_idx,]
      
      cat("Filtered dds object has samples:", ncol(dds_object_filt), "genes:", nrow(dds_object_filt), "\n")
      
      return(dds_object_filt)
      
    } else if (abs_filt == FALSE) {
      # at least 1 fpm in at least x% of samples
      cat(relat_filt * 100, "% of samples:", (relat_filt) * ncol(dds_object), "\n")
      keep_genes_idx <- rowSums(DESeq2::fpm(dds_object, robust = TRUE) > 1) >= relat_filt * ncol(dds_object)
      cat("Number of filtered genes:", sum(keep_genes_idx == FALSE),"\n")
      
      dds_object_filt <- dds_object[keep_genes_idx,]
      
      cat("Filtered dds object has samples:", ncol(dds_object_filt), "genes:", nrow(dds_object_filt), "\n")
      
      return(dds_object_filt)
      
    } else {
      cat("abs_filt should be TRUE or FALSE.\n")
    }
    
  } else {
    cat("Need to specify dds_object; e.g. dds_object=dds \n")
  }
}




#' generatePCA: ADD description
#' 
generatePCA <- function(transf_object=NULL, cond_interest_varPart=NULL, color_variable=NULL, shape_variable=NULL, ntop_genes=500){
  #transf_object_counts <- assay(transf_object)
  #ntop_genes=nrow(transf_object_counts) 
  pcaData <- DESeq2::plotPCA(transf_object, intgroup=cond_interest_varPart, returnData=TRUE, ntop=ntop_genes) 
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  return(
    ggplot(pcaData, aes(PC1, PC2, color=!!sym(color_variable), shape=!!sym(shape_variable))) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + #+ coord_fixed()
    theme_bw() 
  )
}


generateResults <- function(dds_object=NULL, 
                            coef_contrast=NULL,
                            coeff_name=NULL, 
                            cond_numerator=NULL,
                            cond_denominator=NULL, 
                            padj_cutoff = 0.05, 
                            log2FC_cutoff = 0.58,
                            cond_variable=NULL,
                            ensemblAnnot=NULL){
  
  if(!is.null(coef_contrast) & !is.null(coeff_name)){
    stop("Choose only coef_name or coef_contrast")
  }
  
  if(!is.null(coeff_name)){
    print("apeglm shrinkage is used")
    res_table_unshrunken <- DESeq2::results(dds_object, 
                                            name=coeff_name,
                                            parallel = TRUE, alpha = padj_cutoff)
    #renv::install("bioc::apeglm")
    res_table <- DESeq2::lfcShrink(dds_object, 
                                   coef=coeff_name,
                                   res=res_table_unshrunken, type = "apeglm")
  }
  
  if(!is.null(coef_contrast)){
    print("ashr shrinkage is used")
    res_table_unshrunken <- DESeq2::results(dds_object, 
                                            contrast=coef_contrast,
                                            parallel = TRUE, alpha = padj_cutoff)
    #renv::install("bioc::apeglm")
    res_table <- DESeq2::lfcShrink(dds_object, 
                                   contrast=coef_contrast,
                                   res=res_table_unshrunken, type = "ashr")
  }
  
  normalized_counts_AddedMean <- meanExprsPerGroup(dds_object=dds_object,
                                                   cond_numerator=cond_numerator,
                                                   cond_denominator=cond_denominator,
                                                   variable = cond_variable)
  
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("ensembl_id") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "ensembl_id") %>% # adding annotationl entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "ensembl_id") %>% # adding normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(ensembl_id, contains("_symbol"), contains("description"), 
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, contains("stat"), pvalue, padj,
                  gene_biotype,  location,
                  starts_with(paste0(cond_denominator,"-")), starts_with(paste0(cond_numerator,"-"))) %>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    #%>% # everything() - other columns not mentioned; starts_with(paste0(cond_denominator,"_S") may not work if other references names dont follow with _S
    dplyr::rename_with(.fn = ~ gsub("numerator", "", paste0("MeanExpr_", cond_numerator)), .cols = "MeanExpr_numerator") %>% # find better way of renaming!  
    dplyr::rename_with(.fn = ~ gsub("denominator", "", paste0("MeanExpr_", cond_denominator)), .cols = "MeanExpr_denominator")
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = paste0(coeff_name, coef_contrast),
                                        design = paste(as.character(design(dds_object)), collapse=""),
                                        signif_genes = nrow(results_data_annot_signif),
                                        signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
                                        signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff),
                                        cutoffs = paste0("(!is.na(padj) & (padj < ", padj_cutoff,")) & abs(log2FoldChange) > ", log2FC_cutoff))
  results_list<-list(results_signif=results_data_annot_signif, de_details=temp_results_summary_df, results_all=results_data_annot)
  return(results_list)
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
    tibble::rownames_to_column("ensembl_id")
  
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


plotVolcano <- function(dds_results_obj=NULL, genes_of_interest=NULL, plot_title=NULL, log2FC_cutoff=0.58, padj_cutoff=0.05){
  results_data_annot_forPlot <- dds_results_obj %>%
    dplyr::filter(!is.na(pvalue))
  results_data_annot_forPlot$signif_DE <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange > log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results_data_annot_forPlot$signif_DE[results_data_annot_forPlot$log2FoldChange < -log2FC_cutoff & results_data_annot_forPlot$padj < padj_cutoff] <- "DOWN"
  results_data_annot_forPlot$signif_DE <- factor(results_data_annot_forPlot$signif_DE,
                                                 levels = c("NO", "DOWN", "UP"))
  table(results_data_annot_forPlot$signif_DE)
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(padj), col=signif_DE)) +
    geom_point() +
    #gghighlight::gghighlight(signif_DE %in% c("DOWN", "UP")) +
    ggrepel::geom_label_repel(data = . %>% filter(gene_symbol %in% genes_of_interest), 
                              aes(label = gene_symbol),
                              show.legend = FALSE,
                              box.padding = 0.5,
                              segment.color ="black",
                              max.overlaps = Inf
                              ) +
    #geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red", linetype="dashed") +
    #geom_hline(yintercept=-log10(padj_cutoff), col="red", linetype="dashed") + # need to adjust to match padj_cutoff
    #scale_color_manual(values=c(DOWN="navy", UP="firebrick3")) +
    scale_color_manual(values=c(DOWN="navy", UP="firebrick3", NO = "grey")) +
    theme_bw(base_size = 14) +
    labs(x = "log2FC") + 
    ggtitle(plot_title)
  #y = "-log10( p-value )",color = "signif. DE") 
}


generatePCA_repel <- function(transf_object = NULL, cond_interest_varPart = NULL, 
                              color_variable = NULL, shape_variable = NULL, ntop_genes = 500)
{
  pcaData <- DESeq2::plotPCA(transf_object, intgroup = cond_interest_varPart, 
                             returnData = TRUE, ntop = ntop_genes)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, 
                      color = !!sym(color_variable), 
                      shape = !!sym(shape_variable),
                      label = name)) + geom_point(size = 3) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    theme_bw() +
    ggrepel::geom_label_repel(label.padding = 0.1)
}


# ===========================
# ENHANCED UTILITY FUNCTIONS
# ===========================
# Enhanced functions for professional RNA-seq analysis workflows

#' Setup Output Directories
#' @param base_dir Base directory path
#' @param subdirs Vector of subdirectory names to create
setup_output_dirs <- function(base_dir, subdirs = c("plots", "tables", "reports")) {
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
  }
  
  for (subdir in subdirs) {
    subdir_path <- file.path(base_dir, subdir)
    if (!dir.exists(subdir_path)) {
      dir.create(subdir_path, recursive = TRUE)
    }
  }
  
  invisible(TRUE)
}

#' Validate Input Data
#' @param dds DESeq2 object
#' @param annotation_data Annotation data frame
validate_inputs <- function(dds, annotation_data) {
  stopifnot("DESeq2 object is required" = inherits(dds, "DESeqDataSet"))
  stopifnot("Annotation data is required" = is.data.frame(annotation_data))
  stopifnot("ensembl_id column required" = "ensembl_id" %in% colnames(annotation_data))
  stopifnot("gene_symbol column required" = "gene_symbol" %in% colnames(annotation_data))
  
  message("Input validation passed")
  invisible(TRUE)
}

#' Load gene sets for analysis (no caching to avoid repository bloat)
#' @return List of gene sets
load_gene_sets_direct <- function() {
  message("Loading gene sets from MSigDB...")
  gene_sets <- list(
    hallmark = hypeR::msigdb_gsets("Homo sapiens", "H", clean = TRUE),
    kegg = hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean = TRUE),
    reactome = hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean = TRUE),
    gobp = hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean = TRUE),
    gocc = hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean = TRUE),
    gomf = hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean = TRUE),
    oncogenic = hypeR::msigdb_gsets("Homo sapiens", "C6", clean = TRUE),
    wang_hippo = list(genesets = list(wang_hippo_set = c("CCN1", "CCN2", "AMOTL2", "ANKRD1", "IGFBP3", "F3", "FJX1", "NUAK2", "LATS2", "CRIM1", "GADD45A","TGFB2", "PTPN14", "NT5E", "FOXF2", "AXL", "DOCK5", "ASAP1", "RBMS3", "MYOF", "ARHGEF17", "CCDC80")))
  )
  
  return(gene_sets)
}

#' Safely save plots with error handling
#' @param plot_obj ggplot object
#' @param filename Output filename
#' @param ... Additional arguments passed to ggsave
safe_save_plot <- function(plot_obj, filename, ...) {
  tryCatch({
    ggsave(filename = filename, plot = plot_obj, ...)
    message("Saved plot: ", basename(filename))
  }, error = function(e) {
    warning("Failed to save plot ", basename(filename), ": ", e$message)
  })
}

#' Enhanced Volcano Plot with Gene Highlighting
#' @param results DESeq2 results data frame
#' @param genes_highlight Vector of gene symbols to highlight
#' @param title Plot title
#' @param padj_cutoff Adjusted p-value cutoff
#' @param lfc_cutoff Log2 fold change cutoff
create_enhanced_volcano <- function(results, genes_highlight = NULL, title = "Volcano Plot", 
                                  padj_cutoff = 0.05, lfc_cutoff = 1) {
  if (!("gene_symbol" %in% colnames(results))) {
    stop("gene_symbol column required in results")
  }
  
  # Prepare data
  plot_data <- results %>%
    dplyr::mutate(
      significance = dplyr::case_when(
        padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff ~ "Significant",
        TRUE ~ "Not significant"
      ),
      log10_padj = -log10(padj),
      highlight = ifelse(!is.null(genes_highlight) && gene_symbol %in% genes_highlight, "Highlight", "Normal")
    )
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = log2FoldChange, y = log10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = significance, alpha = highlight), size = 0.8) +
    ggplot2::scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey60")) +
    ggplot2::scale_alpha_manual(values = c("Highlight" = 1, "Normal" = 0.6)) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", alpha = 0.5) +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme_publication()
  
  # Add gene labels if specified
  if (!is.null(genes_highlight)) {
    highlight_data <- plot_data %>%
      dplyr::filter(gene_symbol %in% genes_highlight)
    
    if (nrow(highlight_data) > 0) {
      p <- p + 
        ggrepel::geom_text_repel(
          data = highlight_data,
          ggplot2::aes(label = gene_symbol),
          size = 3,
          box.padding = 0.3,
          point.padding = 0.3,
          max.overlaps = 15
        )
    }
  }
  
  return(p)
}

#' Publication-ready Theme
theme_publication <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 11, face = "bold")
    )
}

#' Generate Analysis Summary
#' @param deg_results DESeq results object
#' @param gsea_results Optional GSEA results
generate_analysis_summary <- function(deg_results, gsea_results = NULL) {
  summary_list <- list(
    total_genes_tested = nrow(deg_results$results_all),
    significant_genes = nrow(deg_results$results_signif),
    upregulated = sum(deg_results$results_signif$log2FoldChange > 0),
    downregulated = sum(deg_results$results_signif$log2FoldChange < 0),
    analysis_method = deg_results$de_details$test,
    comparison = paste(deg_results$de_details$numerator, "vs", deg_results$de_details$denominator)
  )
  
  if (!is.null(gsea_results)) {
    summary_list$gsea_analyses <- length(gsea_results)
  }
  
  return(summary_list)
}

#' Enhanced duplicate gene resolution
#' @param deg_results DEG results data frame
#' @param method Method for resolving duplicates
#' @return Data frame with resolved duplicates
resolve_duplicate_genes <- function(deg_results, method = "highest_basemean") {
  initial_count <- nrow(deg_results)
  
  resolved_results <- switch(method,
    "highest_basemean" = deg_results %>%
      group_by(gene_symbol) %>%
      slice_max(baseMean, n = 1, with_ties = FALSE),
    "lowest_pvalue" = deg_results %>%
      group_by(gene_symbol) %>%
      slice_min(pvalue, n = 1, with_ties = FALSE),
    "highest_abs_lfc" = deg_results %>%
      group_by(gene_symbol) %>%
      slice_max(abs(log2FoldChange), n = 1, with_ties = FALSE),
    stop("Invalid method. Use: highest_basemean, lowest_pvalue, or highest_abs_lfc")
  ) %>%
  ungroup()
  
  removed_count <- initial_count - nrow(resolved_results)
  if (removed_count > 0) {
    message("Resolved ", removed_count, " duplicate gene symbols using method: ", method)
  }
  
  return(resolved_results)
}

#' Run GSEA analysis for multiple gene sets
#' @param gene_list Vector of gene symbols
#' @param gene_sets List of gene sets
#' @param background_size Background gene set size
#' @return List of hypeR results
run_gsea_analysis <- function(gene_list, gene_sets, background_size) {
  results <- list()
  
  for (set_name in names(gene_sets)) {
    message("Running GSEA for: ", set_name)
    tryCatch({
      results[[set_name]] <- hypeR::hypeR(
        signature = gene_list,
        genesets = gene_sets[[set_name]],
        test = "hypergeometric",
        background = background_size
      )
    }, error = function(e) {
      warning("GSEA failed for ", set_name, ": ", e$message)
      results[[set_name]] <- NULL
    })
  }
  
  return(results)
}

#' Run fGSEA analysis in batch
#' @param ranked_genes Named numeric vector of ranked genes
#' @param gene_sets List of gene sets
#' @param min_size Minimum gene set size
#' @param max_size Maximum gene set size
#' @param top_n Number of top pathways to return
#' @return List containing fGSEA results and plots
run_fgsea_batch <- function(ranked_genes, gene_sets, min_size = 15, max_size = 500, top_n = 10) {
  results <- list()
  
  for (set_name in names(gene_sets)) {
    message("Running fGSEA for: ", set_name)
    
    tryCatch({
      fgsea_res <- fgsea::fgsea(
        pathways = gene_sets[[set_name]][["genesets"]],
        stats = ranked_genes,
        minSize = min_size,
        maxSize = max_size
      )
      
      # Get top pathways
      top_up <- fgsea_res[ES > 0][head(order(pval), n = top_n), pathway]
      top_down <- fgsea_res[ES < 0][head(order(pval), n = top_n), pathway]
      top_pathways <- c(top_up, rev(top_down))
      
      # Create GSEA table plot
      gsea_plot <- fgsea::plotGseaTable(
        gene_sets[[set_name]][["genesets"]][top_pathways],
        ranked_genes,
        fgsea_res,
        gseaParam = 0.5
      )
      
      # Create dot plot
      dot_plot_data <- fgsea_res %>%
        filter(pathway %in% top_pathways) %>%
        arrange(NES) %>%
        mutate(
          pathway = factor(pathway, levels = pathway),
          padj_log2 = -log2(padj)
        )
      
      dot_plot <- ggplot(dot_plot_data) +
        aes(x = NES, y = pathway, colour = NES, size = padj_log2) +
        geom_point(shape = "circle") +
        scale_color_distiller(palette = "RdBu", direction = -1) +
        scale_size(range = c(4, 10)) +
        theme_publication() +
        labs(
          title = set_name,
          size = "-log2(padj)",
          x = "Normalized Enrichment Score"
        )
      
      results[[set_name]] <- list(
        fgsea_results = fgsea_res,
        gsea_plot = gsea_plot,
        dot_plot = dot_plot,
        top_pathways = top_pathways
      )
      
    }, error = function(e) {
      warning("fGSEA failed for ", set_name, ": ", e$message)
      results[[set_name]] <- NULL
    })
  }
  
  return(results)
}

#' Validate statistical results and provide summary
#' @param results DEG results data frame
#' @param alpha Significance threshold
validate_statistical_results <- function(results, alpha = 0.05) {
  n_tests <- sum(!is.na(results$pvalue))
  bonferroni_threshold <- alpha / n_tests
  
  fdr_significant <- sum(results$padj < alpha, na.rm = TRUE)
  bonf_significant <- sum(results$pvalue < bonferroni_threshold, na.rm = TRUE)
  
  message("=== Statistical Validation ===")
  message("Total tests performed: ", n_tests)
  message("Bonferroni threshold: ", signif(bonferroni_threshold, 3))
  message("FDR significant genes (padj < ", alpha, "): ", fdr_significant)
  message("Bonferroni significant genes: ", bonf_significant)
  message("Median p-value: ", signif(median(results$pvalue, na.rm = TRUE), 3))
  
  return(list(
    n_tests = n_tests,
    bonferroni_threshold = bonferroni_threshold,
    fdr_significant = fdr_significant,
    bonferroni_significant = bonf_significant
  ))
}

message("Enhanced RNA-seq analysis utility functions loaded successfully!")
