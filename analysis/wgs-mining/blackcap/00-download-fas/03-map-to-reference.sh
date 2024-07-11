#!/bin/bash

#SBATCH --job-name=process-sample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --partition=medium
#SBATCH --output=process-sample_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your-email@example.com

module purge
module load SRA-Toolkit/3.0.3-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BWA/0.7.17-GCCcore-11.2.0

set -x
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

sample=$1

# Convert SRA to FASTQ (assuming paired-end data)
fastq-dump --split-files -O $RAW_DATA "$sample"

# Align FASTQ to reference genome and convert directly to BAM
bwa mem $HOST_SYLV_REG "${RAW_DATA}/${sample}_1.fastq" "${RAW_DATA}/${sample}_2.fastq" | samtools view -bS - > "${RAW_DATA}/${sample}.bam"

# Optional: Remove FASTQ files to save space
rm "${RAW_DATA}/${sample}_1.fastq" "${RAW_DATA}/${sample}_2.fastq"


