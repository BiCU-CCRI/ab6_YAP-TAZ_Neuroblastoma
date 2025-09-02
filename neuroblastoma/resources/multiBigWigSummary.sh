#!/bin/bash

multiBigwigSummary BED-file \
--bwfiles ../data/ATACseq/bigWigs/*bigWig \
--BED ../data/ATACseq/bigWigs/CLB_SKN_M.bed \
-out ../data/ATACseq/bigWigs/scores_per_transcript.npz --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript.tab
 
 # multiBigwigSummary BED-file \
 # --bwfiles ../data/ATACseq/bigWigs/*bigWig \
 # --BED ../data/ATACseq/bigWigs/CLB_SKN_M_h3K27Ac.bed \
 # -out ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_MES_SEs_h3k27ac.npz \
 # --outRawCounts ../data/ATACseq/bigWigs/scores_per_transcript_clbm_skn_MES_SEs_h3k27ac.tab
