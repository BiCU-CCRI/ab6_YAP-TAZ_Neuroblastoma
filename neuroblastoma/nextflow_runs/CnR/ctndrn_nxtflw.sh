export username=$(whoami)
NXF_SINGULARITY_CACHEDIR=/home/${username}/singularity_images
export NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR

# design file
input_design=./Cut_and_run_427_263.csv  #for cut and run the format is group, replicate, fastq_1, fastq_2, control
output_dir=./output

# reference files
ref_fasta=./Homo_sapiens_Genome.GRCh38.102.fa
ref_bwa_index=./Human_GRCh38_v102/bowtie2_index
ref_gtf=./Homo_sapiens.GRCh38.102.gtf
ref_bed=./Homo_sapiens.GRCh38.102.bed12
ref_blacklist=./Human_GRCh38_v102/hg38-blacklist.v2.bed

nextflow run nf-core/cutandrun \
    -work-dir /scratch/research/bicu/aleks_b/ab6_soren_cut_and_run/work \
    --normalisation_mode CPM \
    --seacr_norm norm \
    --input ${input_design} \
    --outdir ${output_dir} \
    --genome GRCh38 \
    --bowtie2 ${ref_bwa_index} \
    --gtf ${ref_gtf} \
    --blacklist ${ref_blacklist} \
    --fasta ${ref_fasta} \
    -profile singularity \
    --peakcaller SEACR,MACS2 \
    --normalisation_binsize 1 \
    -c ./biohazard_12c60g.config \
    --dt_calc_all_matrix false 

