#!/bin/bash

# multiBigwigSummary BED-file \
# --bwfiles ../data/ATACseq/bigWigs/*bigWig \
# --BED ../data/ATACseq/bigWigs/CLB_SKN_M.bed \
# -out ../data/ATACseq/bigWigs/scores_per_transcript.npz --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript.tab
 
 # multiBigwigSummary BED-file \
 # --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 # --BED ../data/ATACseq/bigWigs/CLB_SKN_A.bed \
 # -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADR.npz \
 # --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADR.tab
 
 # multiBigwigSummary BED-file \
 # --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 # --BED ../temp_results/BEDs/CLB_SKN_A_u.bed \
 # -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADRu.npz \
 # --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADRu.tab
 # 
 # multiBigwigSummary BED-file \
 # --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 # --BED ../temp_results/BEDs/CLB_SKN_M_u.bed \
 # -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_MESu.npz \
 # --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_MESu.tab
 
 
 multiBigwigSummary BED-file \
 --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 --BED ../temp_results/BEDs/CLB_M_u.bed \
 -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_MESu.npz \
 --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_MESu.tab

 multiBigwigSummary BED-file \
 --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 --BED ../temp_results/BEDs/CLB_A_u.bed \
 -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_ADRu.npz \
 --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_ADRu.tab

 multiBigwigSummary BED-file \
 --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 --BED ../temp_results/BEDs/SKN_M_u.bed \
 -out ../data/ATACseq/bigWigs/scores_per_transcript_skn_MESu.npz \
 --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_skn_MESu.tab

 multiBigwigSummary BED-file \
 --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 --BED ../temp_results/BEDs/SKN_A_u.bed \
 -out ../data/ATACseq/bigWigs/scores_per_transcript_skn_ADRu.npz \
 --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_skn_ADRu.tab