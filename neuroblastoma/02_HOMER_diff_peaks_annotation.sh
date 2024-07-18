#!/bin/bash


homer_dir=/home/rstudio/workspace/neuroblastoma/homer

mkdir ${homer_dir}
wget http://homer.ucsd.edu/homer/configureHomer.pl -O ${homer_dir}/configureHomer.pl
perl ${homer_dir}/configureHomer.pl -install
perl ${homer_dir}/configureHomer.pl -install hg38

PATH=$PATH:/home/rstudio/workspace/neuroblastoma/homer/bin/

perl ${homer_dir}/bin/annotatePeaks.pl \
  /home/rstudio/workspace/neuroblastoma/results/ATAC-seq/bed_diff_ADR.bed \
  hg38 \
  -log \
  -annStats /home/rstudio/workspace/neuroblastoma/results/ATAC-seq/ATAC_genome_states_dataframe_ADR.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
  /home/rstudio/workspace/neuroblastoma/results/ATAC-seq/bed_diff_MES.bed \
  hg38 \
  -log \
  -annStats /home/rstudio/workspace/neuroblastoma/results/ATAC-seq/ATAC_genome_states_dataframe_MES.csv
