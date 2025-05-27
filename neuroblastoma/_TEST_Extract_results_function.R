# test_DDS_extraction <- function (dds_object = NULL, 
#                                  coeff_name = NULL,
#                                  cond_numerator = NULL, 
#                                  cond_denominator = NULL, 
#                                  cond_variable = NULL,
#                                  padj_cutoff = 0.5, 
#                                  log2FC_cutoff = 0.58) 
# {
#   res_table_unshrunken <- DESeq2::results(dds_object,
#                                           name = coeff_name, 
#                                           parallel = TRUE, 
#                                           alpha = padj_cutoff)
#   res_table <- DESeq2::lfcShrink(dds_object, 
#                                  coef = coeff_name, 
#                                  res = res_table_unshrunken, 
#                                  type = "apeglm")
#   row.names(res_table) <- dds_object@rowRanges@elementMetadata@listData[["peak_id"]]
#   normalized_counts_AddedMean <- meanExprsPerGroup(dds_object = dds_object, 
#                                                    cond_numerator = cond_numerator, 
#                                                    cond_denominator = cond_denominator, 
#                                                    variable = cond_variable)
#   normalized_counts_AddedMean$peak_id <- dds_object@rowRanges@elementMetadata@listData[["peak_id"]]
#   ensemblAnnot <- as.data.frame(dds_object@rowRanges@elementMetadata@listData)
#   results_data_annot <- as.data.frame(res_table) %>% tibble::rownames_to_column("peak_id") %>% 
#     dplyr::mutate(FoldChange = ifelse(log2FoldChange < 0, 
#                                       -2^(abs(log2FoldChange)), 2^(log2FoldChange))) %>% 
#     dplyr::left_join(ensemblAnnot, "peak_id") %>% 
#     dplyr::left_join(normalized_counts_AddedMean, "peak_id")  %>% 
#     dplyr::arrange(padj) %>% 
#     dplyr::select(peak_id,
#                   gencode_chr,
#                   gencode_start, 
#                   gencode_end,
#                   is_promoter_3kb,
#                   gene_category_our_RNAseq,
#                   gencode_gene_name,
#                   gencode_characterization,
#                   homer_Nearest.PromoterID, 
#                   homer_Annotation,
#                   homer_Nearest.Ensembl,
#                   homer_Distance.to.TSS,
#                   baseMean.y,
#                   MeanExpr_denominator, 
#                   MeanExpr_numerator,
#                   log2FoldChange,
#                   lfcSE,
#                   FoldChange, 
#                   pvalue,
#                   padj,
#                   starts_with(paste0(cond_denominator, "-")),
#                   starts_with(paste0(cond_numerator, "-"))) %>%
#     dplyr::rename_at(., .vars = "MeanExpr_numerator", 
#                      .funs = funs(gsub("numerator", "", paste0("MeanExpr_", cond_numerator)))) %>% 
#     dplyr::rename_at(., .vars = "MeanExpr_denominator", .funs = funs(gsub("denominator", "", paste0("MeanExpr_", cond_denominator))))
#   results_data_annot_signif <- results_data_annot %>% 
#     dplyr::filter((!is.na(padj) &  (padj < padj_cutoff)) & abs(log2FoldChange) > log2FC_cutoff)
#   temp_results_summary_df <- data.frame(test = coeff_name, 
#                                         design = paste(as.character(design(dds_object)), collapse = ""), 
#                                         signif_genes = nrow(results_data_annot_signif), 
#                                         signif_genes_UP = sum(results_data_annot_signif$log2FoldChange > log2FC_cutoff),
#                                         signif_genes_DOWN = sum(results_data_annot_signif$log2FoldChange < log2FC_cutoff), 
#                                         cutoffs = paste0("(!is.na(padj) &   (padj < ", padj_cutoff, ")) & abs(log2FoldChange) > ", log2FC_cutoff))
#   dds_consens_results_list <- list(results_signif = results_data_annot_signif, 
#                                    de_details = temp_results_summary_df, 
#                                    results_all = results_data_annot)
#   return(dds_consens_results_list)
# }
# 





extract_results_DDS <- function(dds_object = NULL,
                                dds_result = NULL,
                                coeff_name = NULL,
                                cond_numerator = NULL,
                                cond_denominator = NULL,
                                cond_variable = NULL,
                                padj_cutoff = 0.5,
                                log2FC_cutoff = 0.58){
  
  #This is the start of the PR function that I re-purposed.
  if(is.null(dds_result)){
    res_table_unshrunken <- DESeq2::results(dds_object, 
                                            name=coeff_name,   
                                            parallel = TRUE, 
                                            alpha = padj_cutoff)
  }else{
    res_table_unshrunken <- dds_result
  }
  
  res_table <- DESeq2::lfcShrink(dds_object, 
                                 coef=coeff_name,   
                                 res=res_table_unshrunken, 
                                 type = "apeglm")
  row.names(res_table) <-   dds_object@rowRanges@elementMetadata@listData[["PeakId"]] #DATASET   SPECIFIC!
  
  normalized_counts_AddedMean <-   meanExprsPerGroup(dds_object=dds_object,
                                                     cond_numerator=cond_numerator,
                                                     cond_denominator=cond_denominator,
                                                     variable = cond_variable)
  normalized_counts_AddedMean$PeakId <-   dds_object@rowRanges@elementMetadata@listData[["PeakId"]]
  
  ensemblAnnot <- as.data.frame(dds_object@rowRanges@elementMetadata@listData[1:19])
  ######
  #Operations with table, merging with metadata.
  results_data_annot <- as.data.frame(res_table) %>%
    tibble::rownames_to_column("PeakId") %>% 
    dplyr::mutate(FoldChange=ifelse(log2FoldChange < 0, 
                                    -2^(abs(log2FoldChange)), 
                                    2^(log2FoldChange))) %>% #adding FoldChange column
    dplyr::left_join(ensemblAnnot, "PeakId") %>% # adding annotationl   entrez_ids, gene symbols,...
    dplyr::left_join(normalized_counts_AddedMean, "PeakId") %>% # adding   normalized counts and average counts per group
    dplyr::arrange(padj) %>% # ordering based on p.adj
    dplyr::select(PeakId, Chr, Start, End, Strand,
                  Gene.Name, Gene.Type, Annotation, Detailed.Annotation,
                  Distance.to.TSS, Nearest.PromoterID, Entrez.ID,   Nearest.Unigene,
                  baseMean, MeanExpr_denominator, MeanExpr_numerator, log2FoldChange, lfcSE, FoldChange, pvalue, padj,
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
    tibble::rownames_to_column("PeakId")
  
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
