#!/bin/bash
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
neuroblastoma_dir="/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/ab6_YAP-TAZ_Neuroblastoma_restorage/neuroblastoma"

# CnR bam files produced by the nextflow alignment run (bowtie2, marked duplicates)
cnr_input="/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/raw_data/CnR/outdir/02_alignment/bowtie2/target/markdup"

# Output location for deduplicated bams, HMCan calls and LILY super-enhancer output
cnr_lily="${neuroblastoma_dir}/results/CnR/CnR_LILY"
mkdir -p "$cnr_lily"

configFile="${neuroblastoma_dir}/HMCan/configurations/HMCan.config.broad.custom.txt"
fai="${neuroblastoma_dir}/resources/Homo_sapiens_Genome.GRCh38.102.fa.fai"
hmcanBin="${neuroblastoma_dir}/HMCan/src/HMCan"
wigToBigWig="${neuroblastoma_dir}/resources/wigToBigWig"
refseq="${neuroblastoma_dir}/resources/hg38_refseq_noChr.ucsc"
runLily="${neuroblastoma_dir}/LILY/scripts/runLILY.R"

# ---------------------------------------------------------------------------
# Samples
# NAME|SEP, where SEP is the character separating the sample name from the
# mark in the bam file names (some samples use "_", others use "-").
# ---------------------------------------------------------------------------
SAMPLES=(
	"CLB-Ma-A|_"
	"CLB-Ma-M|_"
	"SK-N-SH-A|_"
	"SK-N-SH-M|_"
	"CLB-MA-MES|-"
	"CLB-Ma-ADR|-"
	"SK-N-SH-ADR|-"
	"SK-N-SH-MES|-"
)

# ---------------------------------------------------------------------------
# 1) Remove duplicates (samtools rmdup) for every H3K27ac/IgG replicate bam
#    that exists for each sample.
# ---------------------------------------------------------------------------
# for entry in "${SAMPLES[@]}"; do
# 	name="${entry%%|*}"
# 	sep="${entry##*|}"
# 	for mark in H3K27ac IgG; do
# 		for rep in R1 R2; do
# 			inBam="${cnr_input}/${name}${sep}${mark}_${rep}.target.markdup.sorted.bam"
# 			outBam="${cnr_lily}/${name}_${mark}_${rep}.noDup.bam"
# 			if [ -f "$inBam" ]; then
# 				samtools view -u -q 20 "$inBam" | samtools rmdup -s - "$outBam"
# 			else
# 				echo "NOTE: no ${mark} ${rep} bam for ${name}, skipping (${inBam})"
# 			fi
# 		done
# 	done
# done

# ---------------------------------------------------------------------------
# 2) Run HMCan (H3K27ac vs IgG) + wigToBigWig for every replicate that has an
#    H3K27ac bam. If the matching IgG replicate is missing, fall back to
#    IgG_R1 as background control.
# ---------------------------------------------------------------------------
for entry in "${SAMPLES[@]}"; do
	name="${entry%%|*}"
	for rep in R1 R2; do
		h3k27acBam="${cnr_lily}/${name}_H3K27ac_${rep}.noDup.bam"
		[ -f "$h3k27acBam" ] || continue

		iggBam="${cnr_lily}/${name}_IgG_${rep}.noDup.bam"
		if [ ! -f "$iggBam" ]; then
			iggBam="${cnr_lily}/${name}_IgG_R1.noDup.bam"
			echo "NOTE: no IgG ${rep} bam for ${name}, using IgG_R1 as control"
		fi

		outFolder="${cnr_lily}/${name}_H3K27ac_${rep}/${name}_H3K27ac_${rep}"
		mkdir -p "$(dirname "$outFolder")"
		"$hmcanBin" \
			"$h3k27acBam" \
			"$iggBam" \
			"$configFile" \
			"$outFolder"
		"$wigToBigWig" -clip "${outFolder}.wig" "$fai" "${outFolder}.bw"
	done
done

# ---------------------------------------------------------------------------
# 3) Super-enhancer detection (LILY)
# ---------------------------------------------------------------------------
cat "$runLily" | R --slave --args \
	"${cnr_lily}/" \
	"${cnr_lily}/lily" \
	12500 \
	3000 \
	"$refseq" \
	"$fai"
