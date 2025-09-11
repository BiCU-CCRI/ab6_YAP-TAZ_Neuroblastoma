#!/bin/bash

# multiBigwigSummary BED-file \
# --bwfiles ../data/ATACseq/bigWigs/*bigWig \
# --BED ../data/ATACseq/bigWigs/CLB_SKN_M.bed \
# -out ../data/ATACseq/bigWigs/scores_per_transcript.npz --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript.tab
 
 multiBigwigSummary BED-file \
 --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 --BED ../data/ATACseq/bigWigs/CLB_SKN_A.bed \
 -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADR.npz \
 --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_ADR.tab
