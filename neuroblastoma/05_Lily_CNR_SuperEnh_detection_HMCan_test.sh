configFile="/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/ab6_YAP-TAZ_Neuroblastoma_restorage/neuroblastoma/HMCan/configurations/HMCan.config.broad.custom.txt"
fai="/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/ab6_YAP-TAZ_Neuroblastoma_restorage/neuroblastoma/resources/Homo_sapiens_Genome.GRCh38.102.fa.fai"

cnr_lily="/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/ab6_YAP-TAZ_Neuroblastoma_restorage/neuroblastoma/results/CnR/CnR_LILY"

outFolder="${cnr_lily}/CLB-Ma-A_H3K27ac_R1/CLB-Ma-A_H3K27ac_R1"
mkdir $outFolder
/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/ab6_YAP-TAZ_Neuroblastoma_restorage/neuroblastoma/HMCan/src/HMCan \
	${cnr_lily}/CLB-Ma-A_H3K27ac_R1.noDup.bam \
	${cnr_lily}/CLB-Ma-A_IgG_R1.noDup.bam \
	$configFile \
	$outFolder
./resources/wigToBigWig -clip ${outFolder}.wig $fai ${outFolder}.bw