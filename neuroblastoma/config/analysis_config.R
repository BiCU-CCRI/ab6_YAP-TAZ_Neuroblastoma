# RNA-seq Analysis Configuration File
# Author: Enhanced by Claude Code
# Date: 2025-01-30

# Global analysis configuration
ANALYSIS_CONFIG <- list(
  
  # Project metadata
  project = list(
    name = "YAP_TAZ_Neuroblastoma",
    version = "2.0",
    author = "Aleksandr Bykov",
    description = "YAP/TAZ cooperation with AP-1 in neuroblastoma mesenchymal state"
  ),
  
  # Statistical parameters
  stats = list(
    # Significance thresholds
    padj_cutoff = 0.05,
    pvalue_cutoff = 0.01,
    log2FC_cutoff = 1.0,
    
    # Filtering parameters
    min_counts = 10,
    min_samples = 2,
    
    # Variance explained threshold for PCA
    variance_threshold = 0.6,
    
    # Multiple testing correction
    correction_method = "BH"  # Benjamini-Hochberg
  ),
  
  # Data filtering parameters
  filtering = list(
    # Gene filtering
    abs_filt_samples = 2,
    min_count_threshold = 10,
    max_genes = NULL,  # No limit
    
    # Sample filtering
    exclude_samples = c("STA_NB_6"),  # Known problematic samples
    min_lib_size = 1e6,
    
    # Quality control
    outlier_threshold = 3,  # Standard deviations for outlier detection
    correlation_threshold = 0.8
  ),
  
  # Visualization parameters
  visualization = list(
    # PCA settings
    ntop_genes_pca = 1000,
    pca_components = c(1, 2),
    
    # Plot dimensions
    plot_width = 20,
    plot_height = 20,
    plot_units = "cm",
    plot_dpi = 300,
    
    # Colors
    celltype_colors = c("ADR" = "navy", "MES" = "firebrick3"),
    significance_colors = c("Significant" = "red", "Not significant" = "grey60"),
    
    # GSEA visualization
    top_pathways_show = 20,
    pathway_label_length = 70,
    
    # Heatmap settings
    heatmap_top_genes = 50,
    clustering_distance = "euclidean",
    clustering_method = "complete"
  ),
  
  # Gene set enrichment analysis
  gsea = list(
    # Hypergeometric test settings
    test_type = "hypergeometric",
    background_type = "detected_genes",
    
    # Pathway databases
    databases = list(
      hallmark = TRUE,
      kegg = TRUE,
      reactome = TRUE,
      gobp = TRUE,
      gocc = FALSE,  # Often too many terms
      gomf = FALSE,  # Often too many terms
      oncogenic = TRUE,
      custom = TRUE
    ),
    
    # Filtering
    min_pathway_size = 15,
    max_pathway_size = 500,
    fdr_cutoff = 0.05
  ),
  
  # fGSEA specific parameters
  fgsea = list(
    min_size = 15,
    max_size = 500,
    top_pathways_up = 10,
    top_pathways_down = 10,
    eps = 0,  # For exact p-values
    nproc = 1,  # Number of processors
    gseaParam = 0.5
  ),
  
  # File paths
  paths = list(
    # Input directories
    data_dir = "neuroblastoma/data/RNAseq",
    resources_dir = "neuroblastoma/resources",
    
    # Output directories
    results_dir = "neuroblastoma/results/RNA-seq",
    plots_dir = "neuroblastoma/results/RNA-seq/plots",
    tables_dir = "neuroblastoma/results/RNA-seq/tables",
    reports_dir = "neuroblastoma/results/RNA-seq/reports",
    cache_dir = "neuroblastoma/cache",
    
    # Specific files
    dds_file = "rnaseq_deseq_global_deseq_data_set.rds",
    annotation_file = "rnaseq_deseq_global_annotation_gene.tsv",
    gene_sets_cache = "gene_sets.rds"
  ),
  
  # Genes of interest for highlighting
  genes_of_interest = list(
    # All key genes
    all = c("VIM", "YAP1", "WWTR1", "JUN", "FOSL1", "FOSL2", "PHOX2B", "HAND2", "GATA3"),
    
    # Mesenchymal markers
    mesenchymal = c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1"),
    
    # YAP/TAZ pathway
    yap_taz = c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "LATS1", "LATS2"),
    
    # AP-1 transcription factors
    ap1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3"),
    
    # Adrenergic markers
    adrenergic = c("PHOX2B", "HAND2", "GATA3", "DBH", "TH", "PHOX2A", "ISL1"),
    
    # Neuroblastoma risk genes
    risk_genes = c("MYCN", "ALK", "ATRX", "TP53", "RAS"),
    
    # Cell cycle and proliferation
    proliferation = c("MKI67", "PCNA", "TOP2A", "CCNA2", "CCNB1", "CCND1")
  ),
  
  # Custom gene sets
  custom_gene_sets = list(
    # Wang Hippo pathway genes
    wang_hippo = c(
      "CCN1", "CCN2", "AMOTL2", "ANKRD1", "IGFBP3", "F3", "FJX1", 
      "NUAK2", "LATS2", "CRIM1", "GADD45A", "TGFB2", "PTPN14", 
      "NT5E", "FOXF2", "AXL", "DOCK5", "ASAP1", "RBMS3", "MYOF", 
      "ARHGEF17", "CCDC80"
    ),
    
    # Neuroblastoma cell state genes (from literature)
    nb_mes_signature = c(
      "VIM", "FN1", "CDH2", "SNAI2", "PRRX1", "TWIST1", "ZEB1"
    ),
    
    nb_adr_signature = c(
      "PHOX2B", "HAND2", "GATA3", "DBH", "TH", "PHOX2A", "ISL1", "DLX2"
    )
  ),
  
  # Output formats
  output = list(
    # File formats to generate
    excel = TRUE,
    csv = TRUE,
    tsv = TRUE,
    rds = TRUE,
    
    # Plot formats
    png = TRUE,
    pdf = TRUE,
    svg = FALSE,
    
    # Report formats
    html_report = TRUE,
    pdf_report = FALSE,
    
    # Compression
    compress_results = FALSE
  ),
  
  # Performance settings
  performance = list(
    # Memory management
    max_memory_gb = 16,
    gc_frequency = 10,  # Run garbage collection every N operations
    
    # Parallel processing
    use_parallel = FALSE,
    n_cores = 1,
    
    # Caching
    cache_intermediate = TRUE,
    cache_gene_sets = TRUE
  ),
  
  # Analysis options
  analysis = list(
    # Quality control
    remove_outliers = TRUE,
    check_batch_effects = TRUE,
    
    # Differential expression
    shrink_lfc = TRUE,
    independent_filtering = TRUE,
    
    # Multiple comparisons
    run_pairwise_comparisons = FALSE,
    
    # Advanced analyses
    run_wgcna = FALSE,
    run_pathway_analysis = TRUE,
    run_upstream_regulators = FALSE
  )
)

# Validation function for configuration
validate_config <- function(config) {
  errors <- character(0)
  
  # Check required sections
  required_sections <- c("stats", "filtering", "visualization", "gsea", "paths")
  missing_sections <- setdiff(required_sections, names(config))
  if (length(missing_sections) > 0) {
    errors <- c(errors, paste("Missing config sections:", paste(missing_sections, collapse = ", ")))
  }
  
  # Validate statistical parameters
  if (config$stats$padj_cutoff < 0 || config$stats$padj_cutoff > 1) {
    errors <- c(errors, "padj_cutoff must be between 0 and 1")
  }
  
  if (config$stats$log2FC_cutoff < 0) {
    errors <- c(errors, "log2FC_cutoff must be positive")
  }
  
  # Validate visualization parameters
  if (config$visualization$ntop_genes_pca < 100) {
    errors <- c(errors, "ntop_genes_pca should be >= 100 for meaningful PCA")
  }
  
  # Return validation results
  if (length(errors) > 0) {
    stop("Configuration validation failed:\n", paste(errors, collapse = "\n"))
  } else {
    message("Configuration validation passed")
    return(TRUE)
  }
}

# Function to load and validate configuration
load_analysis_config <- function(config_file = NULL) {
  if (!is.null(config_file) && file.exists(config_file)) {
    source(config_file)
  }
  
  # Validate configuration
  validate_config(ANALYSIS_CONFIG)
  
  # Convert relative paths to absolute paths
  base_path <- here::here()
  for (path_name in names(ANALYSIS_CONFIG$paths)) {
    if (!startsWith(ANALYSIS_CONFIG$paths[[path_name]], "/")) {
      ANALYSIS_CONFIG$paths[[path_name]] <- file.path(base_path, ANALYSIS_CONFIG$paths[[path_name]])
    }
  }
  
  return(ANALYSIS_CONFIG)
}

# Export configuration
if (exists("ANALYSIS_CONFIG")) {
  message("Analysis configuration loaded successfully")
  message("Version: ", ANALYSIS_CONFIG$project$version)
  message("Author: ", ANALYSIS_CONFIG$project$author)
}