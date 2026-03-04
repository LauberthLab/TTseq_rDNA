#!/bin/bash
#SBATCH -A p32170
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 60
#SBATCH -t 04:00:00
#SBATCH --mem=50gb
#SBATCH --job-name="TTseq_pipeline_run"
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err


# run_pipeline.sh — TT-seq pipeline (single SLURM job)
#
# Usage:
#   sbatch run_pipeline.sh                       # Full pipeline
#   sbatch run_pipeline.sh -d                    # Dry-run
#   sbatch run_pipeline.sh -m samples.tsv        # Custom metadata
#   sbatch run_pipeline.sh -s 5 -e 6            # Only steps 5-6
#   sbatch run_pipeline.sh -s 4 -e 4            # Only normalization
#
# Steps:
#   1: Trim (fastp)              5: Downsample (samtools)
#   2: Align hg38+rDNA (STAR)    6: BigWigs — individual (bamCoverage)
#   3: Align dm6 spike-in (STAR) 7: BigWigs — merged reps (bamCoverage)
#   4: Normalize (spike-in)

set -euo pipefail

SCRIPT_DIR="."

# Source config and functions
source "${SCRIPT_DIR}/config.sh"
source "${SCRIPT_DIR}/functions.sh"

# Defaults for step control
START_STEP=${START_STEP:-1}
END_STEP=${END_STEP:-7}

parse_args "$@"
ensure_dirs

# Load modules
module load fastp
module load STAR
module load samtools
module load subread
module load deeptools

# Summary 
NUM_SAMPLES=$(get_sample_count)

log "  TT-seq Pipeline"
log " " 
log "  Metadata:  ${METADATA}"
log "  Samples:   ${NUM_SAMPLES}"
log "  Steps:     ${START_STEP} -> ${END_STEP}"
log "  Dry-run:   ${DRY_RUN}"
log "  Threads:   ${THREADS}"

#Helper: should we run this step?
should_run() {
    local step=$1
    [ "$step" -ge "$START_STEP" ] && [ "$step" -le "$END_STEP" ]
}

# Helper: loop a per-sample function over all samples
run_per_sample() {
    local func=$1
    while IFS= read -r sample; do
        $func "$sample"
    done < <(get_all_samples)
}

# Execute pipeline 

if should_run 1; then
    log "========== STEP 1/7: Trimming =========="
    run_per_sample step_trim
    log "Step 1 complete."
fi

if should_run 2; then
    log "========== STEP 2/7: Align hg38+rDNA =========="
    run_per_sample step_align_hg38
    log "Step 2 complete."
fi

if should_run 3; then
    log "========== STEP 3/7: Align dm6 spike-in =========="
    run_per_sample step_align_spikein
    log "Step 3 complete."
fi

if should_run 4; then
    log "========== STEP 4/7: Spike-in normalization =========="
    step_normalize
    log "Step 4 complete."
fi

if should_run 5; then
    log "========== STEP 5/7: Downsampling =========="
    run_per_sample step_downsample
    log "Step 5 complete."
fi

if should_run 6; then
    log "========== STEP 6/7: Individual bigWigs =========="
    run_per_sample step_bigwig
    log "Step 6 complete."
fi

if should_run 7; then
    log "========== STEP 7/7: Merged bigWigs =========="
    step_merge_bigwigs
    log "Step 7 complete."
fi

log " "
log "  PIPELINE FINISHED"
log "  BigWigs:    ${BW_DIR}"
log "  Norm file:  ${SCRATCH_DIR}/spikein_norm_factors.tsv"
log "  Logs:       ${LOG_DIR}"
log " "
