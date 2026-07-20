#!/bin/bash

# Genome-wide ATAC-seq normalization workflow
# Author: Claude Code Analysis

echo "=== GENOME-WIDE ATAC-seq NORMALIZATION ==="

# Option 1: Use multiBigwigSummary bins for genome-wide coverage
echo "OPTION 1: Genome-wide binned coverage (recommended)"
echo "This samples the entire genome in fixed bins for normalization"

multiBigwigSummary bins \
  --bwfiles ../data/ATACseq/bigWigs/*bigWig \
  --binSize 10000 \
  --numberOfProcessors 4 \
  -out ../data/ATACseq/bigWigs/genome_wide_bins.npz \
  --outRawCounts ../data/ATACseq/bigWigs/genome_wide_bins.tab

echo "Genome-wide binned data saved to: genome_wide_bins.tab"

# Option 2: Use all peaks from peak calling (if available)
echo ""
echo "OPTION 2: All accessible regions (if peak files available)"
echo "This uses all ATAC-seq peaks, not just super enhancers"

# Check if we have consensus peak files
if [ -f "../data/ATACseq/peaks/consensus_peaks.bed" ]; then
    echo "Found consensus peaks file, using for normalization..."
    multiBigwigSummary BED-file \
      --bwfiles ../data/ATACseq/bigWigs/*bigWig \
      --BED ../data/ATACseq/peaks/consensus_peaks.bed \
      -out ../data/ATACseq/bigWigs/all_peaks_scores.npz \
      --outRawCounts ../data/ATACseq/bigWigs/all_peaks_scores.tab
    echo "All peaks data saved to: all_peaks_scores.tab"
else
    echo "No consensus peaks file found at ../data/ATACseq/peaks/consensus_peaks.bed"
    echo "Will proceed with binned approach only"
fi

# Option 3: Calculate scaling factors from bigWig files directly
echo ""
echo "OPTION 3: Direct bigWig scaling factors"
echo "Calculate total coverage from bigWig files"

echo "sample,total_coverage" > ../data/ATACseq/bigWigs/bigwig_scaling_factors.csv
for bigwig in ../data/ATACseq/bigWigs/*.bigWig; do
    filename=$(basename "$bigwig" .bigWig)
    echo "Processing: $filename"
    
    # Use bigWigInfo to get total coverage (if available)
    if command -v bigWigInfo &> /dev/null; then
        total_cov=$(bigWigInfo "$bigwig" | grep "sum of values" | awk '{print $5}')
        echo "$filename,$total_cov" >> ../data/ATACseq/bigWigs/bigwig_scaling_factors.csv
    else
        echo "bigWigInfo not available, skipping direct scaling approach"
        break
    fi
done

echo ""
echo "=== RECOMMENDATIONS ==="
echo "1. Use genome_wide_bins.tab for calculating normalization factors"
echo "2. Apply these factors to your region-specific data (super enhancers)"
echo "3. This ensures unbiased normalization based on genome-wide coverage"
echo ""
echo "Next step: Run genome_wide_normalization.R to calculate proper normalization factors"