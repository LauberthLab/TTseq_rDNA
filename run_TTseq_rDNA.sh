#!/bin/bash
#SBATCH -A <your_alloc>
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 60
#SBATCH -t 10:00:00
#SBATCH --mem=50gb
#SBATCH --job-name="TTseq_pipeline"
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err


# run_pipeline.sh — TT-seq pipeline (single SLURM job)
#
# Usage:
#   sbatch run_pipeline.sh                       # Full pipeline
#   sbatch run_pipeline.sh -d                    # Dry-run
#   sbatch run_pipeline.sh -m samples.tsv        # Custom metadata
#   sbatch run_pipeline.sh -s 5 -e 7            # Only steps 5-7
#   sbatch run_pipeline.sh -s 4 -e 4            # Only normalization
#
# Steps:
#   0: FastQC + FastQ Screen     5: Downsample (samtools)
#   1: Trim (fastp)              6: BigWigs — individual (bamCoverage)
#   2: Align hg38+rDNA (STAR)    7: BigWigs - merged reps (bamCoverage)
#   3: Align dm6 spike-in (STAR) 8: MultiQC
#   4: Normalize (spike-in)


set -euo pipefail

SCRIPT_DIR=./

# Source config and functions
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/functions.sh"

# Defaults for step control
START_STEP=${START_STEP:-0}
END_STEP=${END_STEP:-8}

parse_args "$@"
ensure_dirs

# Load modules
module load fastqc
#module load multiqc
module load fastp
module load STAR
module load samtools
module load subread
#module load deeptools

# Summary
NUM_SAMPLES=$(get_sample_count)

log "============================================"
log "  TT-seq Pipeline"
log "============================================"
log "  Metadata:       ${METADATA}"
log "  Samples:        ${NUM_SAMPLES}"
log "  Steps:          ${START_STEP} -> ${END_STEP}"
log "  Dry-run:        ${DRY_RUN}"
log "  FastQ Screen:   ${RUN_FASTQ_SCREEN}"
log "  Threads:        ${THREADS}"
log "============================================"

# Helpers
should_run() {
    local step=$1
    [ "$step" -ge "$START_STEP" ] && [ "$step" -le "$END_STEP" ]
}

run_per_sample() {
    local func=$1
    while IFS= read -r sample; do
        $func "$sample"
    done < <(get_all_samples)
}

# Execute pipeline
if should_run 0; then
    log "========== STEP 0/8: FastQC + FastQ Screen =========="
    module load fastqc
    step_fastqc
    log "Step 0 complete."

    log "---------- MultiQC (post-QC) ----------"
    module load multiqc
    step_multiqc_qc
    module unload multiqc
    module unload fastqc

fi

if should_run 1; then
    log "========== STEP 1/8: Trimming =========="
    run_per_sample step_trim
    log "Step 1 complete."
fi

if should_run 2; then
    log "========== STEP 2/8: Align hg38+rDNA =========="
    run_per_sample step_align_hg38
    log "Step 2 complete."
fi

if should_run 3; then
    log "========== STEP 3/8: Align dm6 spike-in =========="
    run_per_sample step_align_spikein
    log "Step 3 complete."
fi

if should_run 4; then
    log "========== STEP 4/8: Spike-in normalization =========="
    step_normalize
    log "Step 4 complete."
fi

if should_run 5; then
    log "========== STEP 5/8: Downsampling =========="
    run_per_sample step_downsample
    log "Step 5 complete."
fi

if should_run 6; then
    log "========== STEP 6/8: Individual bigWigs =========="
    module load deeptools
    run_per_sample step_bigwig
    module unload deeptools
    log "Step 6 complete."
fi

if should_run 7; then
    log "========== STEP 7/8: Merged bigWigs =========="
    module load deeptools
    step_merge_bigwigs
    module unload deeptools
    log "Step 7 complete."
fi

if should_run 8; then
    log "========== STEP 8/8: MultiQC (final) =========="
    module load multiqc
    step_multiqc_final
    module unload multiqc
    log "Step 8 complete."
fi

log "============================================"
log "  PIPELINE FINISHED"
log "  QC:              ${QC_DIR}"
log "  MultiQC (QC):    ${MULTIQC_DIR}/multiqc_qc_report.html"
log "  MultiQC (final): ${MULTIQC_DIR}/multiqc_final_report.html"
log "  BigWigs:         ${BW_DIR}"
log "  Norm file:       ${SCRATCH_DIR}/spikein_norm_factors.tsv"
log "  Logs:            ${LOG_DIR}"
log "============================================"
