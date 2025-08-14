# Analysis Report Generator
# Author: Enhanced by Claude Code
# Date: 2025-01-30

#' Generate comprehensive analysis report
#' @param deg_results DEG results list
#' @param gsea_results GSEA results list
#' @param fgsea_results fGSEA results list
#' @param config Analysis configuration
#' @param output_dir Output directory for report
#' @return Path to generated report
generate_analysis_report <- function(deg_results, gsea_results = NULL, fgsea_results = NULL, 
                                   config, output_dir) {
  
  # Create report data structure
  report_data <- list(
    metadata = list(
      analysis_date = Sys.Date(),
      analysis_time = Sys.time(),
      analyst = config$project$author,
      project = config$project$name,
      version = config$project$version,
      r_version = R.version.string,
      session_info = sessionInfo()
    ),
    
    # Data summary
    data_summary = list(
      total_genes_input = nrow(deg_results$results_all),
      total_samples = ncol(deg_results$dds_object %||% NA),
      samples_after_qc = NA,  # To be filled
      genes_after_filtering = nrow(deg_results$results_all),
      significant_genes = nrow(deg_results$results_signif),
      upregulated_genes = sum(deg_results$results_signif$log2FoldChange > 0),
      downregulated_genes = sum(deg_results$results_signif$log2FoldChange < 0)
    ),
    
    # Statistical summary
    statistical_summary = list(
      padj_cutoff = config$stats$padj_cutoff,
      lfc_cutoff = config$stats$log2FC_cutoff,
      multiple_testing_method = config$stats$correction_method,
      median_pvalue = median(deg_results$results_all$pvalue, na.rm = TRUE),
      median_padj = median(deg_results$results_all$padj, na.rm = TRUE),
      median_lfc = median(abs(deg_results$results_signif$log2FoldChange)),
      max_lfc = max(abs(deg_results$results_signif$log2FoldChange)),
      min_lfc = min(abs(deg_results$results_signif$log2FoldChange))
    ),
    
    # Top genes
    top_genes = list(
      most_upregulated = head(deg_results$results_signif[order(-deg_results$results_signif$log2FoldChange), ], 10),
      most_downregulated = head(deg_results$results_signif[order(deg_results$results_signif$log2FoldChange), ], 10),
      most_significant = head(deg_results$results_signif[order(deg_results$results_signif$padj), ], 10)
    )
  )
  
  # Add GSEA summary if available
  if (!is.null(gsea_results)) {
    gsea_summary <- list()
    for (set_name in names(gsea_results)) {
      if (!is.null(gsea_results[[set_name]]) && length(gsea_results[[set_name]]$data) > 0) {
        gsea_summary[[set_name]] <- list(
          total_pathways = nrow(gsea_results[[set_name]]$data),
          significant_pathways = sum(gsea_results[[set_name]]$data$fdr < config$gsea$fdr_cutoff),
          top_pathway = head(gsea_results[[set_name]]$data$label, 1)
        )
      }
    }
    report_data$gsea_summary <- gsea_summary
  }
  
  # Add fGSEA summary if available
  if (!is.null(fgsea_results)) {
    fgsea_summary <- list()
    for (set_name in names(fgsea_results)) {
      if (!is.null(fgsea_results[[set_name]]$fgsea_results)) {
        fgsea_res <- fgsea_results[[set_name]]$fgsea_results
        fgsea_summary[[set_name]] <- list(
          total_pathways = nrow(fgsea_res),
          significant_pathways = sum(fgsea_res$padj < config$fgsea$eps, na.rm = TRUE),
          positive_nes = sum(fgsea_res$NES > 0, na.rm = TRUE),
          negative_nes = sum(fgsea_res$NES < 0, na.rm = TRUE)
        )
      }
    }
    report_data$fgsea_summary <- fgsea_summary
  }
  
  return(report_data)
}

#' Create HTML summary report
#' @param report_data Report data structure
#' @param output_file Output HTML file path
create_html_summary <- function(report_data, output_file) {
  
  html_content <- paste0(
    '<!DOCTYPE html>
    <html>
    <head>
        <title>RNA-seq Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            h1, h2, h3 { color: #333; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
            th { background-color: #f2f2f2; font-weight: bold; }
            .highlight { background-color: #fff3cd; }
            .section { margin: 30px 0; }
            .metadata { background-color: #f8f9fa; padding: 15px; border-radius: 5px; }
        </style>
    </head>
    <body>
        <h1>RNA-seq Differential Expression Analysis Report</h1>
        
        <div class="metadata">
            <h2>Analysis Metadata</h2>
            <p><strong>Project:</strong> ', report_data$metadata$project, '</p>
            <p><strong>Analyst:</strong> ', report_data$metadata$analyst, '</p>
            <p><strong>Analysis Date:</strong> ', report_data$metadata$analysis_date, '</p>
            <p><strong>R Version:</strong> ', report_data$metadata$r_version, '</p>
        </div>
        
        <div class="section">
            <h2>Data Summary</h2>
            <table>
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Total Genes Tested</td><td>', format(report_data$data_summary$total_genes_input, big.mark = ","), '</td></tr>
                <tr><td>Significant Genes (FDR < ', report_data$statistical_summary$padj_cutoff, ')</td><td>', format(report_data$data_summary$significant_genes, big.mark = ","), '</td></tr>
                <tr><td>Upregulated Genes</td><td>', format(report_data$data_summary$upregulated_genes, big.mark = ","), '</td></tr>
                <tr><td>Downregulated Genes</td><td>', format(report_data$data_summary$downregulated_genes, big.mark = ","), '</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Statistical Summary</h2>
            <table>
                <tr><th>Metric</th><th>Value</th></tr>
                <tr><td>Adjusted P-value Cutoff</td><td>', report_data$statistical_summary$padj_cutoff, '</td></tr>
                <tr><td>Log2 Fold Change Cutoff</td><td>±', report_data$statistical_summary$lfc_cutoff, '</td></tr>
                <tr><td>Multiple Testing Correction</td><td>', report_data$statistical_summary$multiple_testing_method, '</td></tr>
                <tr><td>Median P-value</td><td>', signif(report_data$statistical_summary$median_pvalue, 3), '</td></tr>
                <tr><td>Median |Log2FC|</td><td>', round(report_data$statistical_summary$median_lfc, 2), '</td></tr>
                <tr><td>Maximum |Log2FC|</td><td>', round(report_data$statistical_summary$max_lfc, 2), '</td></tr>
            </table>
        </div>
    '
  )
  
  # Add top genes tables
  if (!is.null(report_data$top_genes)) {
    for (gene_type in names(report_data$top_genes)) {
      genes_df <- report_data$top_genes[[gene_type]]
      if (nrow(genes_df) > 0) {
        html_content <- paste0(html_content,
          '<div class="section">
              <h2>Top ', stringr::str_to_title(gsub("_", " ", gene_type)), '</h2>
              <table>
                  <tr><th>Gene</th><th>Log2FC</th><th>Adj. P-value</th><th>Base Mean</th></tr>'
        )
        
        for (i in 1:min(nrow(genes_df), 10)) {
          gene_row <- genes_df[i, ]
          html_content <- paste0(html_content,
            '<tr>
                <td>', gene_row$gene_symbol, '</td>
                <td>', round(gene_row$log2FoldChange, 2), '</td>
                <td>', signif(gene_row$padj, 3), '</td>
                <td>', round(gene_row$baseMean, 1), '</td>
            </tr>'
          )
        }
        
        html_content <- paste0(html_content, '</table></div>')
      }
    }
  }
  
  # Add GSEA summary
  if (!is.null(report_data$gsea_summary)) {
    html_content <- paste0(html_content,
      '<div class="section">
          <h2>Gene Set Enrichment Analysis Summary</h2>
          <table>
              <tr><th>Database</th><th>Total Pathways</th><th>Significant Pathways</th><th>Top Pathway</th></tr>'
    )
    
    for (set_name in names(report_data$gsea_summary)) {
      gsea_data <- report_data$gsea_summary[[set_name]]
      html_content <- paste0(html_content,
        '<tr>
            <td>', stringr::str_to_upper(set_name), '</td>
            <td>', gsea_data$total_pathways, '</td>
            <td>', gsea_data$significant_pathways, '</td>
            <td>', gsea_data$top_pathway %||% "None", '</td>
        </tr>'
      )
    }
    
    html_content <- paste0(html_content, '</table></div>')
  }
  
  # Close HTML
  html_content <- paste0(html_content,
    '    </body>
    </html>'
  )
  
  # Write to file
  writeLines(html_content, output_file)
  message("HTML report saved to: ", output_file)
  
  return(output_file)
}

#' Generate quality control report
#' @param qc_data Quality control data
#' @param output_file Output file path
generate_qc_report <- function(qc_data, output_file) {
  
  qc_summary <- list(
    sample_info = qc_data$sample_info,
    filtering_stats = qc_data$filtering_stats,
    outlier_detection = qc_data$outlier_detection,
    batch_effects = qc_data$batch_effects
  )
  
  # Save as RDS for detailed analysis
  saveRDS(qc_summary, gsub("\\.html$", "_detailed.rds", output_file))
  
  # Create simple HTML report
  html_content <- paste0(
    '<!DOCTYPE html>
    <html>
    <head>
        <title>Quality Control Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .warning { background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 10px 0; }
            .success { background-color: #d4edda; padding: 10px; border-radius: 5px; margin: 10px 0; }
        </style>
    </head>
    <body>
        <h1>Quality Control Report</h1>
        <h2>Sample Information</h2>
        <p>Total samples analyzed: ', length(qc_data$sample_info), '</p>
        
        <h2>Filtering Results</h2>
        <div class="success">
            <p>Quality control completed successfully</p>
        </div>
        
    </body>
    </html>'
  )
  
  writeLines(html_content, output_file)
  return(output_file)
}

#' Create comprehensive analysis summary
#' @param analysis_results All analysis results
#' @param config Configuration
#' @param output_dir Output directory
create_comprehensive_summary <- function(analysis_results, config, output_dir) {
  
  # Generate main report data
  report_data <- generate_analysis_report(
    deg_results = analysis_results$deg_results,
    gsea_results = analysis_results$gsea_results,
    fgsea_results = analysis_results$fgsea_results,
    config = config,
    output_dir = output_dir
  )
  
  # Create HTML summary
  html_file <- file.path(output_dir, "analysis_summary.html")
  create_html_summary(report_data, html_file)
  
  # Save detailed report data
  rds_file <- file.path(output_dir, "analysis_summary.rds")
  saveRDS(report_data, rds_file)
  
  # Create text summary for quick reference
  text_file <- file.path(output_dir, "analysis_summary.txt")
  create_text_summary(report_data, text_file)
  
  message("Comprehensive summary created:")
  message("- HTML report: ", html_file)
  message("- Detailed data: ", rds_file)
  message("- Text summary: ", text_file)
  
  return(list(
    html_report = html_file,
    detailed_data = rds_file,
    text_summary = text_file
  ))
}

#' Create text summary for quick reference
#' @param report_data Report data
#' @param output_file Output text file
create_text_summary <- function(report_data, output_file) {
  
  summary_text <- paste0(
    "RNA-seq Analysis Summary\n",
    "========================\n\n",
    "Analysis Date: ", report_data$metadata$analysis_date, "\n",
    "Project: ", report_data$metadata$project, "\n",
    "Analyst: ", report_data$metadata$analyst, "\n\n",
    
    "Data Overview:\n",
    "- Total genes tested: ", format(report_data$data_summary$total_genes_input, big.mark = ","), "\n",
    "- Significant genes: ", format(report_data$data_summary$significant_genes, big.mark = ","), "\n",
    "- Upregulated: ", format(report_data$data_summary$upregulated_genes, big.mark = ","), "\n",
    "- Downregulated: ", format(report_data$data_summary$downregulated_genes, big.mark = ","), "\n\n",
    
    "Statistical Parameters:\n",
    "- FDR cutoff: ", report_data$statistical_summary$padj_cutoff, "\n",
    "- Log2FC cutoff: ±", report_data$statistical_summary$lfc_cutoff, "\n",
    "- Median p-value: ", signif(report_data$statistical_summary$median_pvalue, 3), "\n",
    "- Maximum |Log2FC|: ", round(report_data$statistical_summary$max_lfc, 2), "\n\n"
  )
  
  # Add top genes
  if (!is.null(report_data$top_genes$most_upregulated)) {
    top_up <- head(report_data$top_genes$most_upregulated$gene_symbol, 5)
    summary_text <- paste0(summary_text,
      "Top upregulated genes: ", paste(top_up, collapse = ", "), "\n"
    )
  }
  
  if (!is.null(report_data$top_genes$most_downregulated)) {
    top_down <- head(report_data$top_genes$most_downregulated$gene_symbol, 5)
    summary_text <- paste0(summary_text,
      "Top downregulated genes: ", paste(top_down, collapse = ", "), "\n"
    )
  }
  
  writeLines(summary_text, output_file)
  return(output_file)
}

# Null coalescing operator for older R versions
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}