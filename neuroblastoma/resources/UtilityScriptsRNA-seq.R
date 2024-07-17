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
                            coeff_name=NULL, 
                            cond_numerator=NULL,
                            cond_denominator=NULL, 
                            padj_cutoff = 0.05, 
                            log2FC_cutoff = 0.58,
                            cond_variable=NULL,
                            ensemblAnnot=NULL){
  
  res_table_unshrunken <- DESeq2::results(dds_object, 
                                          name=coeff_name,
                                          parallel = TRUE, alpha = padj_cutoff)
  #renv::install("bioc::apeglm")
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,
                                 res=res_table_unshrunken, type = "apeglm")
  
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
    dplyr::rename_at(., .vars = "MeanExpr_numerator", .funs = funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% # find better way of renaming!  
    dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs = funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
  
  results_data_annot_signif <- results_data_annot %>%
    dplyr::filter((!is.na(padj) & (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
  
  temp_results_summary_df <- data.frame(test = coeff_name,
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
  
  signif_volcanoPlot <- ggplot(data = results_data_annot_forPlot, aes(x = log2FoldChange, y = -log10(pvalue), col=signif_DE)) +
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
