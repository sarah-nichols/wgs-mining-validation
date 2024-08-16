#!/bin/bash

#SBATCH --job-name=process_samples
#SBATCH --output=process_samples_%A_%a.out
#SBATCH --error=process_samples_%A_%a.err
#SBATCH --clusters=all
#SBATCH --array=1-11
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# Load necessary modules (module names might vary)
module load SRA-Toolkit/3.0.3-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BWA/0.7.17-GCCcore-11.2.0

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

bwa index ${HOST_REF_PATH}/blackcap_ref_reg.fasta

# Assuming accessions.txt is in the current directory and contains one accession number per line
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SAMPLE_ACCESSIONS_SYLV)

# Step 1: Download FASTQ files using fastq-dump
fastq-dump --split-files --gzip --accession ${ACCESSION} --outdir ${BLACKCAP_SAMPLES/raw_data}

# Step 2: Align with bwa and convert to BAM
# Assuming paired-end data. Adjust accordingly if single-end.
bwa mem -t ${SLURM_CPUS_PER_TASK} ${HOST_REF_PATH/blackcap_ref_reg.fasta} \
    ${BLACKCAP_SAMPLES}/raw_data/${ACCESSION}_1.fastq.gz \
    ${BLACKCAP_SAMPLES}/raw_data/${ACCESSION}_2.fastq.gz | \
    samtools view -bS - > ${BLACKCAP_SAMPLES}/raw_data/${ACCESSION}.bam

# Optional: Remove FASTQ files to save space
#rm ${BLACKCAP_SAMPLES}/raw_data/${ACCESSION}_1.fastq.gz ${BLACKCAP_SAMPLES}/raw_data/${ACCESSION}_2.fastq.gz
