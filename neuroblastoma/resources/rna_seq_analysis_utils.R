# RNA-seq Analysis Utility Functions
# Enhanced by Claude Code for improved RNA-seq analysis workflows

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

#' Validate input data for RNA-seq analysis
#' @param dds_object DESeqDataSet object
#' @param annotation_data Gene annotation data frame
validate_inputs <- function(dds_object, annotation_data) {
  # Validate DESeq object
  if (!inherits(dds_object, "DESeqDataSet")) {
    stop("Input must be a DESeqDataSet object")
  }
  
  # Check sample size
  if (ncol(dds_object) < 6) {
    warning("Few samples detected (n=", ncol(dds_object), "), results may be unreliable")
  }
  
  # Validate annotation data
  required_cols <- c("ensembl_id", "gene_symbol")
  missing_cols <- setdiff(required_cols, colnames(annotation_data))
  if (length(missing_cols) > 0) {
    stop("Annotation data missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for missing annotations
  missing_annotations <- sum(is.na(annotation_data$gene_symbol))
  if (missing_annotations > 0) {
    warning(missing_annotations, " genes have missing gene symbols")
  }
  
  message("Input validation passed successfully")
  return(TRUE)
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
  if (!"gene_symbol" %in% colnames(results)) {
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

message("Enhanced RNA-seq analysis utility functions loaded successfully!")

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

#' Publication-ready ggplot theme
#' @return ggplot theme object
theme_publication <- function() {
  theme_bw() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(size = 11, face = "bold")
  )
}

#' Create enhanced volcano plot
#' @param results DEG results data frame
#' @param genes_highlight Vector of genes to highlight
#' @param title Plot title
#' @param padj_cutoff Adjusted p-value cutoff
#' @param lfc_cutoff Log2 fold change cutoff
#' @return ggplot object
create_enhanced_volcano <- function(results, genes_highlight, title, padj_cutoff = 0.05, lfc_cutoff = 1) {
  # Prepare data
  plot_data <- results %>%
    mutate(
      significance = case_when(
        padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff ~ "Significant",
        TRUE ~ "Not significant"
      ),
      highlight = ifelse(gene_symbol %in% genes_highlight, gene_symbol, NA),
      neg_log10_padj = -log10(padj)
    )
  
  # Count significant genes
  sig_up <- sum(plot_data$padj < padj_cutoff & plot_data$log2FoldChange > lfc_cutoff, na.rm = TRUE)
  sig_down <- sum(plot_data$padj < padj_cutoff & plot_data$log2FoldChange < -lfc_cutoff, na.rm = TRUE)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    geom_point(data = filter(plot_data, !is.na(highlight)), 
               color = "red", size = 2, alpha = 0.8) +
    geom_text_repel(aes(label = highlight), 
                    size = 3, max.overlaps = 20, 
                    box.padding = 0.5, point.padding = 0.5) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), 
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(padj_cutoff), 
               linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey60")) +
    theme_publication() +
    labs(
      title = paste0(title, "\n(", sig_up, " upregulated, ", sig_down, " downregulated)"),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme(legend.position = "bottom")
  
  return(p)
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

#' Generate comprehensive analysis summary
#' @param deg_results DEG results list
#' @param gsea_results GSEA results list
#' @return List with summary statistics
generate_analysis_summary <- function(deg_results, gsea_results = NULL) {
  summary_stats <- list(
    total_genes_tested = nrow(deg_results$results_all),
    significant_genes = nrow(deg_results$results_signif),
    upregulated = sum(deg_results$results_signif$log2FoldChange > 0),
    downregulated = sum(deg_results$results_signif$log2FoldChange < 0),
    median_lfc = median(abs(deg_results$results_signif$log2FoldChange)),
    max_lfc = max(abs(deg_results$results_signif$log2FoldChange))
  )
  
  if (!is.null(gsea_results)) {
    # Extract top enriched pathways
    for (set_name in names(gsea_results)) {
      if (!is.null(gsea_results[[set_name]])) {
        hyp_obj <- gsea_results[[set_name]]
        if (length(hyp_obj$data) > 0) {
          top_pathway <- head(hyp_obj$data$label, 1)
          summary_stats[[paste0("top_", set_name)]] <- top_pathway
        }
      }
    }
  }
  
  return(summary_stats)
}

#' Setup output directories safely
#' @param base_dir Base directory path
#' @param subdirs Vector of subdirectory names
setup_output_dirs <- function(base_dir, subdirs = NULL) {
  # Create base directory
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE)
    message("Created output directory: ", base_dir)
  }
  
  # Create subdirectories if specified
  if (!is.null(subdirs)) {
    for (subdir in subdirs) {
      full_path <- file.path(base_dir, subdir)
      if (!dir.exists(full_path)) {
        dir.create(full_path, recursive = TRUE)
        message("Created subdirectory: ", full_path)
      }
    }
  }
}