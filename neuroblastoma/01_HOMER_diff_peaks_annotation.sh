#!/bin/bash

annotatePeaks.pl ./results/20240515/bed_diff_ADR.bed hg38 -logp -annStats ./results/20240515/ATAC_genome_states_dataframe_ADR.csv
annotatePeaks.pl ./results/20240515/bed_diff_MES.bed hg38 -logp -annStats ./results/20240515/ATAC_genome_states_dataframe_MES.csv
