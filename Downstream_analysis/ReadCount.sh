#!/bin/sh
#SBATCH -A p32170
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 60
#SBATCH -t 00:30:00
#SBATCH --mem=50gb
#SBATCH --job-name="fc"


module load subread 

featureCounts \
  -p \
  -T 12 \
  -F SAF \
  -a <path_to counts chrR SAF>\
  -o <path_to counts txt file>\
  <path_to_bams>