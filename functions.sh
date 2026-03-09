#!/bin/bash
# functions.sh — Pipeline step functions
# Per-sample functions take a sample ID argument.
# All-sample functions (normalize, merge, multiqc) iterate internally.

# 
# STEP 0: FastQC + FastQ Screen for a single sample

step_fastqc() {
    log "FASTQC | Running on all FASTQs in ${FASTQ_DIR}"

    run_cmd fastqc \
        ${FASTQ_DIR}/*.fastq.gz \
        -t 30 \
        -o "${QC_DIR}"

    if [ "$RUN_FASTQ_SCREEN" = "true" ]; then
        log "FASTQ_SCREEN | Running on all FASTQs in ${FASTQ_DIR}"

        run_cmd ${FASTQ_SCREEN} \
            ${FASTQ_DIR}/*.fastq.gz \
            --aligner BOWTIE2 \
            --threads 20 \
            --outdir "${QC_DIR}"
    fi
}


# STEP 1: Trim a single sample
step_trim() {
    local sample=$1
    local name=$(get_sample_name "$sample")
    log "TRIM | ${sample} (${name})"

    run_cmd fastp \
        --in1 "${FASTQ_DIR}/${sample}_R1_001.fastq.gz" \
        --in2 "${FASTQ_DIR}/${sample}_R2_001.fastq.gz" \
        --out1 "${TRIM_DIR}/${sample}_clean_R1.fastq.gz" \
        --out2 "${TRIM_DIR}/${sample}_clean_R2.fastq.gz" \
        --thread 12 \
        --detect_adapter_for_pe \
        --html "${TRIM_DIR}/${sample}_fastp_report.html" \
        --json "${TRIM_DIR}/${sample}_fastp_report.json"
}

# STEP 2: Align a single sample to hg38+rDNA

step_align_hg38() {
    local sample=$1
    local name=$(get_sample_name "$sample")
    log "ALIGN hg38+rDNA | ${sample} (${name})"

    run_cmd STAR \
        --runMode alignReads \
        --runThreadN ${THREADS} \
        --genomeDir ${GENOME_DIR} \
        --readFilesIn "${TRIM_DIR}/${sample}_clean_R1.fastq.gz" "${TRIM_DIR}/${sample}_clean_R2.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${sample}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --genomeLoad NoSharedMemory \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignEndsType EndToEnd \
        --outReadsUnmapped Fastx \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic

    run_cmd samtools index -@ ${THREADS} "${ALIGN_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
}

# STEP 3: Align a single sample to dm6 spike-in
step_align_spikein() {
    local sample=$1
    local name=$(get_sample_name "$sample")
    log "ALIGN dm6 spike-in | ${sample} (${name})"

    run_cmd STAR \
        --runMode alignReads \
        --runThreadN ${THREADS} \
        --genomeDir ${SPIKEIN_GENOME_DIR} \
        --readFilesIn "${TRIM_DIR}/${sample}_clean_R1.fastq.gz" "${TRIM_DIR}/${sample}_clean_R2.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${ALIGN_DIR}/${sample}_spikein_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM MD \
        --genomeLoad NoSharedMemory \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignEndsType EndToEnd \
        --outReadsUnmapped Fastx \
        --quantMode TranscriptomeSAM GeneCounts

    run_cmd samtools index -@ ${THREADS} "${ALIGN_DIR}/${sample}_spikein_Aligned.sortedByCoord.out.bam"
}

# STEP 4: Compute spike-in normalization factors (all samples)

step_normalize() {
    log "NORMALIZE | Computing spike-in normalization factors"

    local norm_file="${SCRATCH_DIR}/spikein_norm_factors.tsv"
    declare -A spikein_counts
    local min_count=""

    while IFS= read -r sample; do
        local bam="${ALIGN_DIR}/${sample}_spikein_Aligned.sortedByCoord.out.bam"
        if [ ! -f "$bam" ]; then
            log "ERROR: Spike-in BAM not found: $bam"
            return 1
        fi

        local count
        if [ "$DRY_RUN" = "true" ]; then
            count=1000
            echo "  [DRY-RUN] samtools view -@ ${THREADS} -c -q 255 $bam"
        else
            count=$(samtools view -@ ${THREADS} -c -q 255 "$bam")
        fi

        spikein_counts["$sample"]=$count
        log "  ${sample}: ${count} reads"

        if [ -z "$min_count" ] || [ "$count" -lt "$min_count" ]; then
            min_count=$count
        fi
    done < <(get_all_samples)

    log "  Minimum spike-in count: ${min_count}"

    echo -e "sample_id\tspikein_reads\tnorm_factor" > "$norm_file"
    while IFS= read -r sample; do
        local count=${spikein_counts["$sample"]}
        local frac
        frac=$(echo "scale=5; ${min_count} / ${count}" | bc)
        echo -e "${sample}\t${count}\t${frac}" >> "$norm_file"
        log "  ${sample}: factor=${frac}"
    done < <(get_all_samples)

    log "  Normalization factors -> ${norm_file}"
}

# STEP 5: Downsample a single sample

step_downsample() {
    local sample=$1
    local name=$(get_sample_name "$sample")
    local norm_file="${SCRATCH_DIR}/spikein_norm_factors.tsv"

    if [ ! -f "$norm_file" ]; then
        log "ERROR: Norm factors not found: $norm_file — run step 4 first"
        return 1
    fi

    local frac
    frac=$(awk -F'\t' -v s="$sample" '!/^#/ && $1==s {print $3}' "$norm_file")

    if [ -z "$frac" ]; then
        log "ERROR: No norm factor for ${sample}"
        return 1
    fi

    local input_bam="${ALIGN_DIR}/${sample}_Aligned.sortedByCoord.out.bam"
    local output_bam="${ALIGN_DIR}/${sample}_sub_sorted.bam"

    # If factor >= 1.0 (or very close), just symlink — no downsampling needed
    local skip
    skip=$(awk "BEGIN {print ($frac >= 0.99999) ? 1 : 0}")
    if [ "$skip" -eq 1 ]; then
        log "DOWNSAMPLE | ${sample} (${name}) frac=${frac} — skipping (factor ~1.0)"
        run_cmd cp "$input_bam" "$output_bam"
        run_cmd samtools index -@ ${THREADS} "$output_bam"
        return 0
    fi

    log "DOWNSAMPLE | ${sample} (${name}) frac=${frac}"

    # samtools view -s SEED.FRAC interprets the decimal part as the fraction,
    # so we need to pass e.g. 42.73 for 73%, NOT 42.1 for 100%.
    # Pad fraction to avoid ambiguity: 0.1 -> use as-is (10%), 1.0 -> skipped above.
    local tmp_bam="${ALIGN_DIR}/${sample}_sub.bam"

    run_cmd samtools view -b -@ ${THREADS} \
        -o "$tmp_bam" \
        -s  ${RAND_SEED}.${frac} \
        "$input_bam"

    run_cmd samtools sort -@ ${THREADS} \
        -o "$output_bam" \
        "$tmp_bam"

    run_cmd rm -f "$tmp_bam"
    run_cmd samtools index -@ ${THREADS} "$output_bam"
}
    
# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Generate bigWig for a single sample
# ─────────────────────────────────────────────────────────────────────────────
step_bigwig() {
    local sample=$1
    local name=$(get_sample_name "$sample")
    local input_bam="${ALIGN_DIR}/${sample}_sub_sorted.bam"
    local output_bw="${BW_DIR}/${name}_bs${BINSIZE}.bw"

    if [ ! -f "$input_bam" ]; then
        log "ERROR: Downsampled BAM not found: $input_bam — run step 5 first"
        return 1
    fi

    log "BIGWIG | ${name}"

    run_cmd bamCoverage \
        -b "$input_bam" \
        -o "$output_bw" \
        -p max \
        --normalizeUsing ${NORMALIZE_USING} \
        --binSize ${BINSIZE}
}

# STEP 7: Merge replicates per condition, generate merged bigWigs (all samples)

step_merge_bigwigs() {
    log "MERGE BIGWIGS | Merging replicates per condition"

    while IFS= read -r group; do
        log "  Merge group: ${group}"

        local bam_list=()
        while IFS= read -r sample; do
            local bam="${ALIGN_DIR}/${sample}_sub_sorted.bam"
            if [ ! -f "$bam" ] && [ "$DRY_RUN" != "true" ]; then
                log "ERROR: BAM not found: $bam"
                return 1
            fi
            bam_list+=("$bam")
        done < <(get_samples_in_group "$group")

        local merged_bam="${ALIGN_DIR}/${group}_merged.bam"
        local merged_bw="${BW_DIR}/${group}_merged_bs${BINSIZE}.bw"

        log "  Merging ${#bam_list[@]} BAMs -> ${merged_bam}"
        run_cmd samtools merge -f -@ ${THREADS} "$merged_bam" "${bam_list[@]}"
        run_cmd samtools index -@ ${THREADS} "$merged_bam"

        log "  BigWig -> ${merged_bw}"
        run_cmd bamCoverage \
            -b "$merged_bam" \
            -o "$merged_bw" \
            -p max \
            --normalizeUsing ${NORMALIZE_USING} \
            --binSize ${BINSIZE}

    done < <(get_merge_groups)
}

# MULTIQC — post-QC (runs after step 0, covers FastQC + FastQ Screen only)
step_multiqc_qc() {
    log "MULTIQC (post-QC) | Aggregating QC reports"

    run_cmd multiqc \
        "${QC_DIR}" "${TRIM_DIR}" \
        -o "${MULTIQC_DIR}" \
        -n multiqc_qc_report \
        --force

    log "  Post-QC report -> ${MULTIQC_DIR}/multiqc_qc_report.html"
}

# STEP 8: MultiQC — final (runs at end, covers everything)

step_multiqc_final() {
    log "MULTIQC (final) | Aggregating all reports from ${SCRATCH_DIR}"

    run_cmd multiqc \
        "$SCRATCH_DIR" \
        -o "${MULTIQC_DIR}" \
        -n multiqc_final_report \
        --force

    log "  Final report -> ${MULTIQC_DIR}/multiqc_final_report.html"
}


step_poli_efficiency() {
    local ttseq_bw=$1
    local chip_bw=$2
    local output_bw=$3

    if [ ! -f "$ttseq_bw" ]; then
        log "ERROR: TT-seq bigWig not found: $ttseq_bw"
        return 1
    fi

    if [ ! -f "$chip_bw" ]; then
        log "ERROR: Pol I ChIP bigWig not found: $chip_bw"
        return 1
    fi

    log "POLI EFFICIENCY | $(basename "$ttseq_bw") / $(basename "$chip_bw")"

    run_cmd bigwigCompare \
        --bigwig1 "$ttseq_bw" \
        --bigwig2 "$chip_bw" \
        --binSize ${BINSIZE} \
        -p max \
        -o "$output_bw"

    log "  -> ${output_bw}"
}
