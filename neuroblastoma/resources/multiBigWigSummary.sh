#!/bin/bash
set -euo pipefail

# deepTools multiBigwigSummary driver (BED-file mode).
#
# deepTools is installed in a dedicated conda env at /home/abykov/envs/deeptools
# (mamba create -p /home/abykov/envs/deeptools -c conda-forge -c bioconda deeptools).
#
# Two modes:
#   1. No arguments  -> default MES-vs-ADRN cell-line comparisons: scores all
#      bigWigs in BW_DIR over the CLB/SKN unique-DAR region BEDs, writing four
#      score tables into OUT_DIR.
#   2. Custom single run:
#        multiBigWigSummary.sh <BED> <OUT_TAB> [<bigWig> ...]
#      Scores the given bigWigs (default: all in BW_DIR) over <BED>, writing
#      <OUT_TAB> (--outRawCounts) plus a sibling .npz. Used by
#      06_ATAC-seq_K975_SEs.Rmd to score MES super-enhancer accessibility.
#
# NOTE: BED chromosome names must match the bigWigs. The nf-core/atacseq
# bigWigs here use bare names (1, 2, ... X), so use the NON-"chr" BEDs.

# --- Paths (edit here) -------------------------------------------------------
DEEPTOOLS_BIN=/home/abykov/envs/deeptools/bin
BW_DIR=/nobackup/lab_ccri_bicu/internal/abykov/projects/ab6_soeren_neuroblastoma/adaptation_for_cemm/raw_data/ATAC-seq/outdir/bwa/merged_library/bigwig
BED_DIR=/home/abykov/neuroblastoma/temp_results/BEDs
OUT_DIR=/home/abykov/neuroblastoma/results/ATAC-seq/ATAC-seq_MES_vs_ADR/multiBigwigSummary

NPROC=8   # number of processors for multiBigwigSummary
# -----------------------------------------------------------------------------

export PATH="${DEEPTOOLS_BIN}:${PATH}"
shopt -s nullglob

run_summary () {
    # run_summary <bed_path> <out_tab> <bigWig> [<bigWig> ...]
    local bed="$1" out_tab="$2"; shift 2
    local out_npz="${out_tab%.tab}.npz"
    mkdir -p "$(dirname "$out_tab")"
    echo ">>> multiBigwigSummary: $(basename "$bed") -> ${out_tab}  ($# bigWigs)"
    multiBigwigSummary BED-file \
        --bwfiles "$@" \
        --BED "$bed" \
        --numberOfProcessors "${NPROC}" \
        -out "${out_npz}" \
        --outRawCounts "${out_tab}"
}

# --- Custom single-run mode --------------------------------------------------
if [[ $# -ge 1 ]]; then
    if [[ $# -lt 2 ]]; then
        echo "Usage: $(basename "$0") <BED> <OUT_TAB> [<bigWig> ...]" >&2
        exit 2
    fi
    bed="$1"; out_tab="$2"; shift 2
    if [[ $# -eq 0 ]]; then
        set -- "${BW_DIR}"/*.bigWig   # default to all bigWigs
    fi
    run_summary "$bed" "$out_tab" "$@"
    echo "Done. Output: ${out_tab}"
    exit 0
fi

# --- Default mode: MES-vs-ADRN cell-line comparisons over all bigWigs --------
mkdir -p "${OUT_DIR}"
BWFILES=( "${BW_DIR}"/*.bigWig )
# NOTE: this includes the drug-treatment (CM/SH/k975) and OE samples as well as
# the MES/ADRN cell lines. Restrict BWFILES if you only want the cell lines.

run_summary "${BED_DIR}/CLB_M_u.bed" "${OUT_DIR}/scores_per_transcript_clbm_MESu.tab" "${BWFILES[@]}"
run_summary "${BED_DIR}/CLB_A_u.bed" "${OUT_DIR}/scores_per_transcript_clbm_ADRu.tab" "${BWFILES[@]}"
run_summary "${BED_DIR}/SKN_M_u.bed" "${OUT_DIR}/scores_per_transcript_skn_MESu.tab"  "${BWFILES[@]}"
run_summary "${BED_DIR}/SKN_A_u.bed" "${OUT_DIR}/scores_per_transcript_skn_ADRu.tab"  "${BWFILES[@]}"

run_summary "${BED_DIR}/CLB_SKN_M_u.bed" "${OUT_DIR}/scores_per_transcript_CLB_SKN_M_u.tab"  "${BWFILES[@]}"

echo "Done. Outputs written to ${OUT_DIR}"
