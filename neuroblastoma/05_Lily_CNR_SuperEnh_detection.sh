#!/bin/bash

# IMPORTANT - Here, the location of the cut and run bam files that got produced during the nextflow run should be specified
cnr_input="~/workspace/neuroblastoma/nextflow_runs/CnR/output"
cnr_lily="~/workspace/neuroblastoma/nextflow_runs/CnR/CnR_LILY/"

samtools view -u -q 20 ${cnr_input}/CLB-Ma-A_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-A_IgG_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/CLB-Ma-A_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-A_IgG_R2.noDup.bam
samtools view -u -q 20 ${cnr_input}/CLB-Ma-A_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-A_H3K27ac_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/CLB-Ma-A_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-A_H3K27ac_R2.noDup.bam

samtools view -u -q 20 ${cnr_input}/CLB-Ma-M_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-M_IgG_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/CLB-Ma-M_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-M_H3K27ac_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/CLB-Ma-M_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/CLB-Ma-M_H3K27ac_R2.noDup.bam

samtools view -u -q 20 ${cnr_input}/SK-N-SH-A_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-A_IgG_R1.target.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-A_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-A_IgG_R2.target.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-A_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-A_H3K27ac_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-A_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-A_H3K27ac_R2.noDup.bam

samtools view -u -q 20 ${cnr_input}/SK-N-SH-M_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-M_IgG_R1.target.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-M_IgG_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-M_IgG_R2.target.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-M_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-M_H3K27ac_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/SK-N-SH-M_H3K27ac_R2.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/SK-N-SH-M_H3K27ac_R2.noDup.bam

samtools view -u -q 20 ${cnr_input}/Ex251-STA-NB-10-MES-100k_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/Ex251-STA-NB-10-MES-100k_IgG_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/Ex251-STA-NB-10-MES-200k_IgG_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/Ex251-STA-NB-10-MES-200k_IgG_R1.noDup.bam

samtools view -u -q 20 ${cnr_input}/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.noDup.bam
samtools view -u -q 20 ${cnr_input}/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.target.markdup.sorted.bam | samtools rmdup -s - ${cnr_lily}/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.noDup.bam


configFile="~/workspace/neuroblastoma/HMCan/configurations/HMCan.config.broad.custom.txt"
fai="~/workspace/neuroblastoma/resources/Homo_sapiens_Genome.GRCh38.102.fa.fai"

outFolder="${cnr_lily}/CLB-Ma-A_H3K27ac_R1/CLB-Ma-A_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/CLB-Ma-A_H3K27ac_R1.noDup.bam \
	${cnr_lily}/CLB-Ma-A_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

#
outFolder="${cnr_lily}/CLB-Ma-A_H3K27ac_R2/CLB-Ma-A_H3K27ac_R2"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/CLB-Ma-A_H3K27ac_R2.noDup.bam \
	${cnr_lily}/CLB-Ma-A_IgG_R2.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/CLB-Ma-M_H3K27ac_R1/CLB-Ma-M_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/CLB-Ma-M_H3K27ac_R1.noDup.bam \
	${cnr_lily}/CLB-Ma-M_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/CLB-Ma-M_H3K27ac_R2/CLB-Ma-M_H3K27ac_R2"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/CLB-Ma-M_H3K27ac_R2.noDup.bam \
	${cnr_lily}/CLB-Ma-M_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/SK-N-SH-A_H3K27ac_R1/SK-N-SH-A_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/SK-N-SH-A_H3K27ac_R1.noDup.bam \
	${cnr_lily}/SK-N-SH-A_IgG_R1.target.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/SK-N-SH-A_H3K27ac_R2/SK-N-SH-A_H3K27ac_R2"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/SK-N-SH-A_H3K27ac_R2.noDup.bam \
	${cnr_lily}/SK-N-SH-A_IgG_R2.target.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/SK-N-SH-M_H3K27ac_R1/SK-N-SH-M_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/SK-N-SH-M_H3K27ac_R1.noDup.bam \
	${cnr_lily}/SK-N-SH-M_IgG_R1.target.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/SK-N-SH-M_H3K27ac_R2/SK-N-SH-M_H3K27ac_R2"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/SK-N-SH-M_H3K27ac_R2.noDup.bam \
	${cnr_lily}/SK-N-SH-M_IgG_R2.target.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/Ex251-STA-NB-10-MES-100k_H3K27ac_R1/Ex251-STA-NB-10-MES-100k_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/Ex251-STA-NB-10-MES-100k_H3K27ac_R1.noDup.bam \
	${cnr_lily}/Ex251-STA-NB-10-MES-100k_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw

outFolder="${cnr_lily}/Ex251-STA-NB-10-MES-200k_H3K27ac_R1/Ex251-STA-NB-10-MES-200k_H3K27ac_R1"
mkdir $outFolder
./HMCan/src/HMCan \
	${cnr_lily}/Ex251-STA-NB-10-MES-200k_H3K27ac_R1.noDup.bam \
	${cnr_lily}/Ex251-STA-NB-10-MES-200k_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw


fai="~/workspace/neuroblastoma/resources/Homo_sapiens_Genome.GRCh38.102.fa.fai"
cat ~/workspace/neuroblastoma/LILY/scripts/runLILY.R | R --slave --args \
        ${cnr_lily}/ \
        ${cnr_lily}/lily \
        12500 \
        3000 \
        ~/workspace/neuroblastoma/resources/hg38_refseq.ucsc \
        $fai