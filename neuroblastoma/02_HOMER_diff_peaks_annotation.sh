#!/bin/bash

# Set HOMER directory
homer_dir=~/neuroblastoma/homer

# HOMER's scripts call helper tools (bed2pos.pl, assignGenomeAnnotation, ...)
# by bare name, so homer/bin must be on PATH
export PATH="${homer_dir}/bin:${PATH}"

perl ${homer_dir}/bin/annotatePeaks.pl \
  ~/neuroblastoma/results/ATAC-seq/ATAC-seq_MES_vs_ADR/bed_diff_ADR.bed \
  hg38 \
  -log \
  -annStats ~/neuroblastoma/results/ATAC-seq/ATAC-seq_MES_vs_ADR/ATAC_genome_states_dataframe_ADR.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
  ~/neuroblastoma/results/ATAC-seq/ATAC-seq_MES_vs_ADR/bed_diff_MES.bed \
  hg38 \
  -log \
  -annStats ~/neuroblastoma/results/ATAC-seq/ATAC-seq_MES_vs_ADR/ATAC_genome_states_dataframe_MES.csv
