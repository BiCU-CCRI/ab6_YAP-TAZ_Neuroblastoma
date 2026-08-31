#!/bin/bash
# Set HOMER directory
homer_dir=~/neuroblastoma/homer

# HOMER's scripts call helper tools (bed2pos.pl, assignGenomeAnnotation, ...)
# by bare name, so homer/bin must be on PATH
export PATH="${homer_dir}/bin:${PATH}"

perl ${homer_dir}/bin/annotatePeaks.pl \
~/neuroblastoma/temp_results/BEDs/CLB_SKN_A_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_A_homer_annot.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
 ~/neuroblastoma/temp_results/BEDs/CLB_SKN_M_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_M_homer_annot.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
 ~/neuroblastoma/temp_results/BEDs/CLB_SKN_AM_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_AM_homer_annot.csv

# annotate the _u - super unique SEs
perl ${homer_dir}/bin/annotatePeaks.pl \
~/neuroblastoma/temp_results/BEDs/CLB_SKN_A_u_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_A_u_homer_annot.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
 ~/neuroblastoma/temp_results/BEDs/CLB_SKN_M_u_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_M_u_homer_annot.csv

perl ${homer_dir}/bin/annotatePeaks.pl \
 ~/neuroblastoma/temp_results/BEDs/CLB_SKN_AM_u_chr.bed hg38 > ~/neuroblastoma/temp_results/BEDs/CLB_SKN_AM_u_homer_annot.csv