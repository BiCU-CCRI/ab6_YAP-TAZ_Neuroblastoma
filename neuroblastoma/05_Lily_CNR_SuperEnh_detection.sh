#!/bin/bash

# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-A_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-A_IgG_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-A_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-A_IgG_R2.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-A_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-A_H3K27ac_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-A_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-A_H3K27ac_R2.noDup.bam
#
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-M_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-M_IgG_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-M_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-M_H3K27ac_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/CLB-Ma-M_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/CLB-Ma-M_H3K27ac_R2.noDup.bam
#
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-A_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-A_IgG_R1.target.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-A_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-A_IgG_R2.target.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-A_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-A_H3K27ac_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-A_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-A_H3K27ac_R2.noDup.bam
#
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-M_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-M_IgG_R1.target.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-M_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-M_IgG_R2.target.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-M_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-M_H3K27ac_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/SK-N-SH-M_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/SK-N-SH-M_H3K27ac_R2.noDup.bam
#
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/Ex251-STA-NB-10-MES-100k_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/Ex251-STA-NB-10-MES-100k_IgG_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/Ex251-STA-NB-10-MES-200k_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/Ex251-STA-NB-10-MES-200k_IgG_R1.noDup.bam
#
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.noDup.bam
# samtools view -u -q 20 ./data_soren/Cut_and_run_ppln_427_263/output_CMT_normalisation/02_alignment/markdup/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ./CnR_LILY/test/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.noDup.bam

# configFile="./hmcan/configurations/HMCan.config.broad.custom.txt"
fai=./resources/Homo_sapiens_Genome.GRCh38.102.fa.fai
#
outFolder="./CnR_LILY/20240604/CLB-Ma-A_H3K27ac_R1/CLB-Ma-A_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/CLB-Ma-A_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/CLB-Ma-A_IgG_R1.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

#
outFolder="./CnR_LILY/20240604/CLB-Ma-A_H3K27ac_R2/CLB-Ma-A_H3K27ac_R2"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/CLB-Ma-A_H3K27ac_R2.noDup.bam \
# 	./CnR_LILY/test/CLB-Ma-A_IgG_R2.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/CLB-Ma-M_H3K27ac_R1/CLB-Ma-M_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/CLB-Ma-M_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/CLB-Ma-M_IgG_R1.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/CLB-Ma-M_H3K27ac_R2/CLB-Ma-M_H3K27ac_R2"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/CLB-Ma-M_H3K27ac_R2.noDup.bam \
# 	./CnR_LILY/test/CLB-Ma-M_IgG_R1.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/SK-N-SH-A_H3K27ac_R1/SK-N-SH-A_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/SK-N-SH-A_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/SK-N-SH-A_IgG_R1.target.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/SK-N-SH-A_H3K27ac_R2/SK-N-SH-A_H3K27ac_R2"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/SK-N-SH-A_H3K27ac_R2.noDup.bam \
# 	./CnR_LILY/test/SK-N-SH-A_IgG_R2.target.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/SK-N-SH-M_H3K27ac_R1/SK-N-SH-M_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/SK-N-SH-M_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/SK-N-SH-M_IgG_R1.target.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/SK-N-SH-M_H3K27ac_R2/SK-N-SH-M_H3K27ac_R2"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/SK-N-SH-M_H3K27ac_R2.noDup.bam \
# 	./CnR_LILY/test/SK-N-SH-M_IgG_R2.target.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/Ex251-STA-NB-10-MES-100k_H3K27ac_R1/Ex251-STA-NB-10-MES-100k_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/Ex251-STA-NB-10-MES-100k_IgG_R1.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#
outFolder="./CnR_LILY/20240604/Ex251-STA-NB-10-MES-200k_H3K27ac_R1/Ex251-STA-NB-10-MES-200k_H3K27ac_R1"
# mkdir $outFolder
# ./HMCan/src/HMCan \
# 	./CnR_LILY/test/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.noDup.bam \
# 	./CnR_LILY/test/Ex251-STA-NB-10-MES-200k_IgG_R1.noDup.bam \
# 	$configFile \
# 	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw
#

# fai=/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Public/references/Human_GRCh38_v102/Homo_sapiens_Genome.GRCh38.102.fa.fai
# cat ./LILY/scripts/runLILY.R | R --slave --args \
# 	/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/CnR_LILY/test/nhCan/ \
# 	/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/CnR_LILY/test/nhCan/lily \
# 	12500 \
# 	3000 \
# 	/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Internal/aleksandr/projects/ab6_20230508_soren_YAP_TAZ/strohmenger/neuroblastoma/resources/hg38_refseq.ucsc \
# 	/home/aleksandr_b/bioinf_isilon/core_bioinformatics_unit/Public/references/Human_GRCh38_v102/Homo_sapiens_Genome.GRCh38.102.fa.fai
