#!/bin/bash

bw1="./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/SH_A_ATAC.bigWig"
bw2="./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/SH_M_ATAC.bigWig"
bed1="./tmp/tmp_bed_cluster_1.bed"
bed2="./tmp/tmp_bed_cluster_2.bed"

#bigwigAverage -b \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/SH_A_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/NB8_A_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/CLBM_A_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/NB10_A_ATAC.bigWig \
#	-o \
#	./tmp/test_merged_averaged_A_ATAC.bigWig

#bigwigAverage -b \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/SH_M_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/NB8_M_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/CLBM_M_ATAC.bigWig \
#	./data_soren/ATAC_BSA_0729_KB_NB_plasticity/atacseq_hub/NB10_M_ATAC.bigWig \
#	-o \
#	./tmp/test_merged_averaged_M_ATAC.bigWig

computeMatrix scale-regions -S \
	./tmp/test_merged_averaged_A_ATAC.bigWig \
	./tmp/test_merged_averaged_M_ATAC.bigWig \
	-R $bed1 $bed2 \
	-o ./tmp/test_compute_matrix_merged.gz -m 5000
plotHeatmap -m ./tmp/test_compute_matrix_merged.gz -o ./tmp/test_compute_matrix.png
