#!/bin/sh
#SBATCH -A p31270
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 70
#SBATCH -t 2:00:00
#SBATCH --mem=50gb
#SBATCH --job-name="PolI_TTseq"



source config.sh
source functions.sh
module load deeptools


step_poli_efficiency \
    ${BW_DIR}/TT_DMSO_merged_bs1.bw \
    /path/to/chip/DMSO_PolI_merged.bw \
    ${BW_DIR}/DMSO_PolI_efficiency.bw

step_poli_efficiency \
    ${BW_DIR}/TT_ISD_merged_bs1.bw \
    /path/to/chip/ISD_PolI_merged.bw \
    ${BW_DIR}/ISD_PolI_efficiency.bw