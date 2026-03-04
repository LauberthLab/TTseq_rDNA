# TT-seq rDNA Pipeline

A modular, spike-in normalized TT-seq pipeline for mapping nascent transcription across the human ribosomal DNA locus. Designed for SLURM-based HPC clusters (tested on Northwestern Quest).

---

## Overview

This pipeline processes paired-end TT-seq data through trimming, alignment to a custom hg38+rDNA reference genome, *Drosophila* spike-in normalization, and bigWig generation - both per-replicate and merged across conditions. It is built around a metadata-driven design: a single TSV file defines all samples, conditions, and replicate groupings.

### Pipeline Steps

| Step | Description | Tool | Mode |
|------|-------------|------|------|
| 0 | FastQC + FastQ Screen (optional) | FastQC, FastQ Screen | per-sample |
| 1 | Adapter trimming | fastp | per-sample |
| 2 | Align to hg38+rDNA | STAR | per-sample |
| 3 | Align to dm6 spike-in | STAR | per-sample |
| 4 | Compute spike-in normalization factors | samtools + bc | all samples |
| 5 | Downsample BAMs by norm factor | samtools | per-sample |
| 6 | Generate individual replicate bigWigs | deeptools bamCoverage | per-sample |
| 7 | Merge replicates per condition, generate merged bigWigs | samtools merge + bamCoverage | all samples |
| 8 | MultiQC (final) | MultiQC | all samples |

MultiQC runs twice: once automatically after Step 0 (covering FastQC/FastQ Screen results for a quick quality check) and again at Step 8 (aggregating all reports across the full pipeline).

### Dependency Graph

```
Step 0 (fastqc) → MultiQC (post-QC report)
Step 1 (trim adapters and barcodes)
  ├──→ Step 2 (align hg38) ──────────┐
  └──→ Step 3 (align dm6) → Step 4 (normalize) ──→ Step 5 (downsample) ──→ Step 6 (bigwigs) ──┐
                                       Step 2 ──────────────────────────↗         └──→ Step 7 (merge bigwigs) ──→ Step 8 (multiqc final)
```

---

## Prerequisites

### Reference Genome: hg38 + rDNA

Standard hg38 lacks a complete ribosomal DNA unit. This pipeline uses a custom reference that appends the full human rDNA repeat (GenBank U13369.1, ~44 kb) to the hg38 assembly.

**Obtaining the genome files:**

Pre-built FASTA and GTF files optimized for rDNA mapping are available from:

> [https://github.com/vikramparalkar/rDNA-Mapping-Genomes](https://github.com/vikramparalkar/rDNA-Mapping-Genomes)

The default ones are already added in the `Annotation` directory

This repository provides `hg38+rDNA` FASTA files and corresponding GTF annotations that include proper gene features across the 47S precursor (5'ETS, 18S, ITS1, 5.8S, ITS2, 28S, 3'ETS) as well as the intergenic spacer (IGS).

### Building the STAR Index

Once you have the combined FASTA and GTF, build a STAR genome index:

```bash
STAR \
  --runMode genomeGenerate \
  --runThreadN 16 \
  --genomeDir /path/to/STAR_genome+rDNA \
  --genomeFastaFiles /path/to/genome+rDNA fasta  \
  --sjdbGTFfile /path/to/genome+rDNA GTF \
  --sjdbOverhang <median read length-1>
```
Set `--sjdbOverhang` to your median read length minus 1 (e.g., 99 for 100 bp reads).

The default ones on Quest: 
FASTA files: /projects/b1042/LauberthLab/Genome/Human_hg38-rDNA_genome_v1.0
GTF files: /projects/b1042/LauberthLab/Genome/20220806_hg38_refseq_chrR_KY962518.1_35500cut_renamed.gtf


Similarly, build a STAR index for the *Drosophila* dm6 spike-in genome using the standard dm6 FASTA and GTF from Ensembl or UCSC.

### Pre-built indices on Quest (Northwestern)

If you are on Quest, pre-built indices are already available:

| Genome | Path |
|--------|------|
| hg38+rDNA | `/projects/b1042/LauberthLab/Genome/STAR_hg38+rDNA` |
| dm6 (spike-in) | `/projects/b1042/LauberthLab/Genome/STAR_index_dm6` |

These paths are set as defaults in `config.sh`.

### Required Modules (Quest)

```
fastqc (0.12.0)
multiqc (1.14)
fastp (0.23.4)
STAR (2.7.5a)
samtools (1.6)
subread (2.0.3)
deeptools (3.5.6)
FastQ-Screen(0.15.3)
```

FastQ Screen is optional and uses a local install at `/projects/b1042/LauberthLab/FastQ-Screen-0.15.3/fastq_screen`. Enable it with the `-f` flag.

All are available via `module load` on Quest and are loaded automatically by the pipeline.

---

## Project Structure

```
TTseq_rDNA_pipeline/
├── config.sh              # Paths, parameters, metadata helpers
├── functions.sh           # Step functions (trim, align, normalize, etc.)
├── run_TTseq_rDNA.sh      # Main SLURM script (single-job sequential mode)
├── metadata.tsv           # Sample metadata (edit this for your experiment)
├── DE_analysis.R           # Differential expression analysis (DESeq2)
├── PolI_efficiency_BigWig.sh   # Pol I elongation efficiency bigWig generation
├── Pol1_ChIP_comparison.py     # Pol I ChIP-seq comparison analysis
├── PolI_ChIP_comparison.ipynb  # Jupyter notebook for ChIP comparison
├── Cleavage_analysis.py        # rRNA cleavage site analysis
├── ReadCount.sh                # Read counting utility
└── logs/                       # Auto-created log directory
```

---

## Quick Start

### 1. Prepare your metadata

Create or edit `metadata.tsv` — a tab-separated file with four columns:

```
#sample_id	condition	replicate	merge_group
21161R-97-01_S96_L002	DMSO	rep1	TT_DMSO
21161R-97-02_S97_L002	DMSO	rep2	TT_DMSO
21161R-97-03_S98_L002	DMSO	rep3	TT_DMSO
21161R-97-04_S99_L002	TPT	rep1	TT_ISD
21161R-97-05_S100_L002	TPT	rep2	TT_ISD
21161R-97-06_S101_L002	TPT	rep3	TT_ISD
```

| Column | Description |
|--------|-------------|
| `sample_id` | FASTQ prefix (expects `{sample_id}_R1_001.fastq.gz` and `_R2_001.fastq.gz`) |
| `condition` | Experimental condition label |
| `replicate` | Replicate identifier (rep1, rep2, ...) |
| `merge_group` | Samples sharing a merge_group are combined in Step 7 for merged bigWigs |

### 2. Update paths in config.sh

Edit `config.sh` to point to your data:

```bash
SCRATCH_DIR="/projects/b1042/LauberthLab/TOP1_nucleolus/TTseq_rDNA_pipeline/Data"
FASTQ_DIR="${SCRATCH_DIR}/FASTQ"         # Where your .fastq.gz files live
```

The FASTQ directory should contain files named as:
```
{sample_id}_R1_001.fastq.gz
{sample_id}_R2_001.fastq.gz
```

### 3. Run the full pipeline

```bash
sbatch run_TTseq_rDNA.sh
```

### 4. Monitor

```bash
squeue -u $USER
tail -f logs/pipeline_*.out
```

---

## Usage Options

```bash
# Full pipeline
sbatch run_TTseq_rDNA.sh

# Dry-run (prints commands without executing)
sbatch run_TTseq_rDNA.sh -d

# Custom metadata file
sbatch run_TTseq_rDNA.sh -m /path/to/my_samples.tsv

# Enable FastQ Screen (off by default)
sbatch run_TTseq_rDNA.sh -f

# Run only specific steps
sbatch run_TTseq_rDNA.sh -s 5 -e 7      # Downsample + bigWigs only
sbatch run_TTseq_rDNA.sh -s 4 -e 4      # Recompute normalization only
sbatch run_TTseq_rDNA.sh -s 0 -e 0      # QC only (FastQC + MultiQC)
sbatch run_TTseq_rDNA.sh -s 6 -e 8      # Regenerate bigWigs + final MultiQC
```

| Flag | Description | Default |
|------|-------------|---------|
| `-m` | Path to metadata TSV | `./metadata.tsv` |
| `-d` | Dry-run mode | off |
| `-f` | Enable FastQ Screen | off |
| `-s` | Start from step N | 0 |
| `-e` | End after step N | 8 |
| `-h` | Print help | — |

---

## Pipeline Details

### Quality Control (Step 0)

Runs FastQC on raw FASTQ files (both R1 and R2) for each sample, with reports written to `QC/`. Optionally runs FastQ Screen (enabled with `-f`) to check for contamination from other species. A post-QC MultiQC report (`multiqc_qc_report.html`) is generated automatically after this step so you can assess read quality before committing to the full pipeline.

### Trimming (Step 1)

Uses fastp with automatic adapter detection for paired-end reads. Generates per-sample HTML and JSON QC reports in the `FASTQ_TRIM/` directory.

### Alignment (Steps 2–3)

STAR is run in two-pass mode for the primary hg38+rDNA alignment (Step 2) and single-pass for the dm6 spike-in (Step 3). Key STAR parameters:

- `--outFilterMultimapNmax 10` — allows up to 10 multimapping positions (important for rDNA, which has ~300–400 tandem copies), balances for the fact that the reads still align to a single annotation
- `--outFilterMismatchNoverReadLmax 0.04` — strict mismatch filter
- `--alignEndsType EndToEnd` — no soft-clipping

### Spike-in Normalization (Step 4)

Counts uniquely mapped reads (MAPQ = 255) in each dm6 spike-in BAM. Computes a normalization fraction as `min_count / sample_count` for each sample. The sample with the fewest spike-in reads gets a factor of 1.0; all others are downsampled proportionally. Factors are saved to `spikein_norm_factors.tsv`.

### Downsampling (Step 5)

Uses `samtools view -s` with the computed normalization fraction and a fixed random seed (42) for reproducibility. The output is a sorted, indexed BAM.

### BigWig Generation (Steps 6–7)

Uses deeptools `bamCoverage` with RPKM normalization and 1 bp bin size. Step 6 generates per-replicate bigWigs from downsampled BAMs. Step 7 merges all replicates within each `merge_group` (defined in metadata) using `samtools merge`, then generates a merged bigWig per condition.

### MultiQC (Step 8)

Runs MultiQC across the entire `SCRATCH_DIR`, aggregating reports from FastQC, FastQ Screen, fastp, STAR, and deeptools into a single final report (`multiqc_final_report.html`). This complements the post-QC report generated after Step 0, giving a comprehensive view of the entire pipeline run.

---

## Output Files

After a successful run, key outputs are organized as:

```
Data/
├── QC/
│   ├── {sample}_R1_001_fastqc.html
│   ├── {sample}_R2_001_fastqc.html
│   └── {sample}_R1_001_screen.html          # if -f enabled
├── MULTIQC/
│   ├── multiqc_qc_report.html               # post-QC (after Step 0)
│   └── multiqc_final_report.html            # final (after Step 8)
├── FASTQ_TRIM/
│   ├── {sample}_clean_R1.fastq.gz
│   ├── {sample}_clean_R2.fastq.gz
│   └── {sample}_fastp_report.html
├── BAMS/
│   ├── {sample}_Aligned.sortedByCoord.out.bam       # hg38+rDNA alignment
│   ├── {sample}_spikein_Aligned.sortedByCoord.out.bam  # dm6 alignment
│   ├── {sample}_sub_sorted.bam                       # downsampled BAM
│   ├── {sample}_ReadsPerGene.out.tab                  # STAR gene counts
│   └── {merge_group}_merged.bam                       # merged BAM
├── BIGWIGS/
│   ├── {merge_group}_{rep}_bs1.bw                     # per-replicate
│   └── {merge_group}_merged_bs1.bw                    # merged
└── spikein_norm_factors.tsv                           # normalization table
```

---

## Downstream Analysis Scripts

The following scripts are included for downstream analyses (some are under active development):

| Script | Description |
|--------|-------------|
| `DE_analysis.R` | Differential expression analysis using DESeq2 on STAR gene counts |
| `PolI_efficiency_BigWig.sh` | Generate Pol I elongation efficiency tracks from TT-seq data |
| `Pol1_ChIP_comparison.py` | Compare TT-seq signal with Pol I ChIP-seq occupancy |
| `PolI_ChIP_comparison.ipynb` | Interactive Jupyter notebook for Pol I ChIP comparison |
| `Cleavage_analysis.py` | Analyze rRNA cleavage patterns from nascent RNA data |
| `ReadCount.sh` | Utility script for counting reads across rDNA features |

---

## Configuration Reference

All tunable parameters are in `config.sh`:

| Variable | Default | Description |
|----------|---------|-------------|
| `SCRATCH_DIR` | `.../TTseq_rDNA_pipeline/Data` | Root data directory |
| `GENOME_DIR` | `.../STAR_hg38+rDNA` | STAR index for hg38+rDNA |
| `SPIKEIN_GENOME_DIR` | `.../STAR_index_dm6` | STAR index for dm6 |
| `FASTQ_SCREEN` | `.../FastQ-Screen-0.15.3/fastq_screen` | Path to FastQ Screen executable |
| `THREADS` | 40 | Threads for STAR, samtools, bamCoverage |
| `RAND_SEED` | 42 | Random seed for downsampling reproducibility |
| `BINSIZE` | 1 | bigWig bin size (bp) |
| `NORMALIZE_USING` | RPKM | deeptools normalization method |
| `RUN_FASTQ_SCREEN` | false | Enable FastQ Screen (set via `-f` flag) |

---

## Troubleshooting

**Pipeline fails at Step 4 (normalize):** Ensure all dm6 spike-in alignments completed (Step 3). Check that spike-in BAMs exist in `BAMS/` with the `_spikein_` prefix.

**BigWigs look empty:** Verify that `spikein_norm_factors.tsv` has reasonable values (all between 0 and 1). A factor near 0 means that sample had vastly more spike-in reads than others — check for contamination or library prep issues.

**STAR runs out of memory:** The hg38+rDNA index requires ~32 GB RAM. Ensure your SLURM allocation has `--mem=50gb` or higher.

**Module not found:** On Quest, ensure you're on a compute node (not a login node) and that the module names haven't changed. Run `module avail <tool>` to check.

---

## References

- Paralkar, V. et al. — rDNA Mapping Genomes: [https://github.com/vikramparalkar/rDNA-Mapping-Genomes](https://github.com/vikramparalkar/rDNA-Mapping-Genomes)
- Schwalb, B. et al. (2016). TT-seq maps the human transient transcriptome. *Science*, 352(6290), 1225–1228.
- Dobin, A. et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21.
- Ramírez, F. et al. (2016). deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research*, 44(W1), W160–W165.
