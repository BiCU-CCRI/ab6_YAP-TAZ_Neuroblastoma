# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview
This is a bioinformatics research project analyzing YAP/TAZ cooperation with AP-1 transcription factors in driving the mesenchymal cell state in neuroblastoma. The project processes ATAC-seq, RNA-seq, and Cut&Run data to investigate cellular plasticity and therapy resistance mechanisms.

## Development Environment Setup
This project uses a Docker-based R development environment with RStudio Server.

### Initial Setup Commands
```bash
# 1. Build the Docker environment
./00_create_environment.sh

# 2. Launch RStudio Server
./run_rstudio.sh

# 3. Connect to RStudio at localhost:48907 (or IP_ADDRESS:48907)
#    Username: rstudio
#    Password: test0

# 4. Open workspace.Rproj in RStudio

# 5. Setup the project environment
source("neuroblastoma/00_setup_the_project.R")
```

## Key Analysis Scripts and Execution Order
The analysis should be run in this specific order:

### Core Analysis Scripts
1. `neuroblastoma/00_setup_the_project.R` - Environment setup, package installation, directory creation
2. `neuroblastoma/01_RNA_seq_analysis.Rmd` - Primary RNA-seq analysis
3. `neuroblastoma/01_RNA_seq_OverExpression_analysis.Rmd` - RNA-seq overexpression analysis
4. `neuroblastoma/01_RNA-seq_YAP_TAZ_inhibition.Rmd` - YAP/TAZ inhibition analysis
5. `neuroblastoma/01_RNA_seq_WWTR1_OE_analysis.Rmd` - WWTR1 overexpression analysis
6. `neuroblastoma/02_ATAC_seq_analysis.R` - ATAC-seq data processing
7. `neuroblastoma/03_ATAC_seq_analysis.Rmd` - ATAC-seq analysis notebook
8. `neuroblastoma/02_CnR_analysis.Rmd` - Cut&Run analysis
9. `neuroblastoma/03_CNR_analysis.R` - CNR processing scripts

### Supporting Analysis Scripts  
- `neuroblastoma/06_RNA-seq_OE_analysis.R` - Additional overexpression analysis
- `neuroblastoma/06_RNA-seq_yap_taz_inhibition.R` - YAP/TAZ inhibition processing
- `neuroblastoma/07_soeren_heatmaps_NB_plasticity_RNAseq.R` - Heatmap generation
- `neuroblastoma/08_Supplementary_visualisations.R` - Supplementary figures

## Project Architecture

### Data Organization
- `neuroblastoma/data/` - Input data files organized by assay type
  - `ATACseq/` - ATAC-seq data and processed objects
  - `RNAseq/` - RNA-seq DESeq objects and annotations  
  - `CnR/` - Cut&Run peak data
  - `other/` - Clinical data and gene sets
  - `temp_data/` - Temporary analysis files

### Results Structure
- `neuroblastoma/results/` - Analysis outputs organized by assay
  - `RNA-seq/` - Differential expression results
  - `ATAC-seq/` - Chromatin accessibility results and plots
  - `CnR/` - Cut&Run peak enrichment analyses

### Third-party Tools (Submodules)
- `neuroblastoma/HMCan/` - Histone modification detection tool (C++)
- `neuroblastoma/LILY/` - Super-enhancer detection (R scripts)  
- `neuroblastoma/homer/` - Motif analysis and peak annotation toolkit

## R Environment Management
- Uses `renv` for package management and reproducibility
- R version 4.2.0 with Bioconductor 3.16
- Package dependencies defined in `renv.lock`
- Packages automatically restored via `renv::restore()` in setup script

## Data Processing Workflow
1. **RNA-seq**: DESeq2-based differential expression analysis of MES vs ADR phenotypes
2. **ATAC-seq**: ChromVAR analysis for chromatin accessibility and transcription factor activity
3. **Cut&Run**: Peak calling and motif enrichment for YAP/TAZ and AP-1 binding sites
4. **Integration**: Cross-assay analysis to identify core regulatory circuits

## Key Utility Scripts
- `neuroblastoma/resources/utilityScripts.R` - General utility functions
- `neuroblastoma/resources/UtilityScriptsRNA-seq.R` - RNA-seq specific utilities  
- `neuroblastoma/00_load_data_and_modules.sh` - Data loading automation

## Build Commands for Submodules
```bash
# Compile HMCan (if needed)
cd neuroblastoma/HMCan/src && make

# HOMER is installed automatically via setup script
```

## Common Development Tasks
- Knit R Markdown files to generate analysis reports with embedded plots
- Use `renv::snapshot()` to update package lockfile after adding dependencies
- Results are automatically saved to appropriate subdirectories in `neuroblastoma/results/`