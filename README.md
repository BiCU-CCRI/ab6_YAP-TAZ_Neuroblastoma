# YAP/TAZ cooperate with AP-1 to drive the mesenchymal cell state in neuroblastoma

## Overview
This code was used for processing ATAC-seq, RNA-seq and Cut&Run data obtained in this study.

## Abstract
Cellular plasticity is a major driver of therapy resistance in cancer. For the pediatric solid tumor neuroblastoma, two distinct tumor cell phenotypes have been identified: adrenergic (ADR) and mesenchymal (MES) cells. Of those, the MES cells pose a rare subtype resistant to common treatment options. Here, we developed an image-based compound screen to interrogate vulnerabilities specific to the MES phenotype at single-cell resolution in heterogeneous neuroblastoma cultures. Among the top hits, we identified inhibition of the transcriptional co-activators YAP/TAZ as a specific target in MES cells. Based on their chromatin binding, we show that YAP/TAZ cooperate with AP-1 transcription factors to regulate MES gene expression by forming a core-regulatory circuitry. We further provide evidence of a subset of tumor cells in patients co-expressing these factors. Collectively, we demonstrate that YAP/TAZ and AP-1 drive neuroblastoma tumor cell plasticity and present a novel therapeutic vulnerability with the potential to overcome treatment resistance.

## System requirements 
### Software requirements 
**OS Requirements**
This code was developed and tested on Linux Debian 8 (jessie)
This code requires Docker version 18.06.3-ce, build d7080c1 to launch R-Server or local R Studio with R v4.2.0
R dependencies are listed in renv.lock file and managed by renv

## Instrucitons of use
Code split into individual blocks, each for specific part of the analysis. And supposed to be run in this particular 

## HOW TO USE DEVELOPMENT
1) Clone the repository, including submodules with 
`git clone --recurse-submodules https://github.com/BiCU-CCRI/ab6_YAP-TAZ_Neuroblastoma`
2) Run `00_create_environment.sh` to create a docker image from a docker file.
3) Run the `run_rstudio.sh`
4) Connect to the Rstudio using "ipaddress.of.your.machine:48907" and loging in to the R-Server with "rstudio" Username and "test0" Password
5) Open a project workspace.Rproj
6) Setup the project (folders, packages, etc. using `00_setup_the_project.R`)
This will do several things: 
  - Reconstruct the R environment using renv framework.
  - Create folders that are used for storing results
7) Run the analysises in the following order:
  - 01_RNA_seq_analysis.Rmd 
  - 01_RNA_seq_OverExpression_analysis.Rmd



# TODOs
01_RNA_seq_analysis.Rmd
- Is the export of mes_adrn_genes to a list actually necessary??


## 3rd party tools used:
wigToBigWig https://www.encodeproject.org/software/wigtobigwig/
HMCan used as a submodule
LILY used as a submodule



