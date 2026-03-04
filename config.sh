#!/bin/bash
###############################################################################
# config.sh — Shared configuration for TT-seq pipeline
###############################################################################

# ── Paths ────────────────────────────────────────────────────────────────────
SCRATCH_DIR="/projects/b1042/LauberthLab/TOP1_nucleolus/TTseq_rDNA_pipeline/Data"
FASTQ_DIR="${SCRATCH_DIR}/FASTQ"
TRIM_DIR="${SCRATCH_DIR}/FASTQ_TRIM"
ALIGN_DIR="${SCRATCH_DIR}/BAMS"
BW_DIR="${SCRATCH_DIR}/BIGWIGS"
GENOME_DIR="/projects/b1042/LauberthLab/Genome/STAR_hg38+rDNA"
SPIKEIN_GENOME_DIR="/projects/b1042/LauberthLab/Genome/STAR_index_dm6"

# ── Parameters ───────────────────────────────────────────────────────────────
THREADS=40
RAND_SEED=42
BINSIZE=1
NORMALIZE_USING="RPKM"

# ── Metadata (override with -m flag) ─────────────────────────────────────────
METADATA="${METADATA:-$(dirname "${BASH_SOURCE[0]}")/metadata.tsv}"

# ── Dry-run mode ─────────────────────────────────────────────────────────────
DRY_RUN="${DRY_RUN:-false}"

run_cmd() {
    if [ "$DRY_RUN" = "true" ]; then
        echo "  [DRY-RUN] $*"
    else
        echo "  [RUN] $*"
        "$@"
    fi
}

# ── Logging ──────────────────────────────────────────────────────────────────
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${PIPELINE_DIR}/logs"
mkdir -p "${LOG_DIR}" 2>/dev/null

log() {
    local ts
    ts=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${ts}] $*" | tee -a "${LOG_DIR}/pipeline_$(date +%Y%m%d).log"
}

# ── Metadata parsing ─────────────────────────────────────────────────────────

get_all_samples() {
    awk -F'\t' '!/^#/ && NF>=1 {print $1}' "$METADATA"
}

get_sample_count() {
    get_all_samples | wc -l
}

get_field() {
    local sample_id=$1 field=$2
    awk -F'\t' -v s="$sample_id" -v f="$field" '!/^#/ && $1==s {print $f}' "$METADATA"
}

get_sample_name() {
    local sample_id=$1
    echo "$(get_field "$sample_id" 4)_$(get_field "$sample_id" 3)"
}

get_merge_groups() {
    awk -F'\t' '!/^#/ && NF>=4 {print $4}' "$METADATA" | sort -u
}

get_samples_in_group() {
    local group=$1
    awk -F'\t' -v g="$group" '!/^#/ && $4==g {print $1}' "$METADATA"
}

# ── Argument parsing ─────────────────────────────────────────────────────────

parse_args() {
    while getopts "m:ds:e:h" opt; do
        case $opt in
            m) METADATA="$OPTARG" ;;
            d) DRY_RUN="true" ;;
            s) START_STEP=$OPTARG ;;
            e) END_STEP=$OPTARG ;;
            h)
                echo "Usage: $0 [-m metadata.tsv] [-d] [-s start_step] [-e end_step] [-h]"
                echo "  -m  Path to metadata file (default: ./metadata.tsv)"
                echo "  -d  Dry-run mode (print commands without executing)"
                echo "  -s  Start from step N (default: 1)"
                echo "  -e  End after step N (default: 7)"
                echo "  -h  Show this help"
                echo ""
                echo "Steps: 1=trim, 2=align_hg38, 3=align_spikein,"
                echo "       4=normalize, 5=downsample, 6=bigwigs, 7=merge_bigwigs"
                exit 0
                ;;
            *) echo "Unknown option" >&2; exit 1 ;;
        esac
    done

    if [ ! -f "$METADATA" ]; then
        echo "ERROR: Metadata file not found: $METADATA" >&2
        exit 1
    fi
}

ensure_dirs() {
    mkdir -p "$TRIM_DIR" "$ALIGN_DIR" "$BW_DIR" 2>/dev/null
}
