#!/bin/bash

multiBigwigSummary BED-file \
 --bwfiles ./*bigWig \
 --BED ./CLB_SKN_M.bed \
 -out ./scores_per_transcript.npz --outRawCounts ./scores_per_transcript.tab