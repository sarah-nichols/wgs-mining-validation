#!/bin/bash

#SBATCH --job-name=process_samples
#SBATCH --output=process_samples_%A_%a.out
#SBATCH --error=process_samples_%A_%a.err
#SBATCH --array=1-11
#SBATCH --mem=10G
#SBATCH --cpus-per-task=4

# Load necessary modules (module names might vary)
module load SRA-Toolkit/3.0.3-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BWA/0.7.17-GCCcore-11.2.0

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

# Assuming accessions.txt is in the current directory and contains one accession number per line
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_ACCESSIONS_SYLV)

# Step 1: Download FASTQ files using fastq-dump
fastq-dump --split-files --gzip --accession ${ACCESSION} --outdir ${RAW_DATA}

# Step 2: Align with bwa and convert to BAM
# Assuming paired-end data. Adjust accordingly if single-end.
bwa mem -t ${SLURM_CPUS_PER_TASK} ${$HOST_SYLV_REG} \
    ${RAW_DATA}/${ACCESSION}_1.fastq.gz \
    ${RAW_DATA}/${ACCESSION}_2.fastq.gz | \
    samtools view -bS - > ${RAW_DATA}/${ACCESSION}.bam

# Optional: Remove FASTQ files to save space
rm ${RAW_DATA}/${ACCESSION}_1.fastq.gz ${RAW_DATA}/${ACCESSION}_2.fastq.gz
