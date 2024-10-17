#!/bin/bash

#SBATCH --job-name=process_samples
#SBATCH --output=process_samples_%A_%a.out
#SBATCH --error=process_samples_%A_%a.err
#SBATCH --clusters=all
#SBATCH --array=1-11
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

# Load necessary modules
module load SRA-Toolkit/3.0.3-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BWA/0.7.17-GCCcore-11.2.0
module load fastp/0.23.4-GCC-12.2.0

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <parameters_file> [--mail-user=email@example.com]"
    exit 1
fi

# Assign input arguments to variables
PARAMETERS_FILE=$1

# Optionally accept email as a command-line argument (if provided, it will be used in the job submission)
EMAIL=${2:-}

# Set the email if provided
if [ -n "$EMAIL" ]; then
    sbatch --mail-user="$EMAIL"
fi

# Load parameters from the file
declare -A PARAMS
while IFS="=" read -r key value; do
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$PARAMETERS_FILE"

# Retrieve the required parameters from the file
ACCESSION_LIST=${PARAMS["accession_list"]}
BLACKCAP_SAMPLES=${PARAMS["blackcap_samples"]}
HOST_REF_PATH=${PARAMS["host_index_path"]}
REPORT_DIR=${PARAMS["report_dir"]}

# Validate that all required parameters are provided
if [ -z "$ACCESSION_LIST" ] || [ -z "$BLACKCAP_SAMPLES" ] || [ -z "$HOST_REF_PATH" ] || [ -z "$REPORT_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Get the current accession from the accession list based on the SLURM_ARRAY_TASK_ID
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}{p;q}" "$ACCESSION_LIST" | tr -d '\r' | tr -d '\n')

# Debugging output: print the retrieved accession
echo "Retrieved accession: '$ACCESSION'"

# Ensure the ACCESSION is not empty
if [ -z "$ACCESSION" ]; then
    echo "Error: Accession number could not be retrieved for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Define paths for raw data and filtered data
RAW_DATA_DIR="${BLACKCAP_SAMPLES}/raw_data"
FILTERED_DATA_DIR="${BLACKCAP_SAMPLES}/filtered_data"

# Create necessary directories if they don't exist
#mkdir -p "$RAW_DATA_DIR"
mkdir -p "$FILTERED_DATA_DIR"
mkdir -p "$REPORT_DIR"

# Step 1: Download FASTQ files using fastq-dump
#echo "Downloading FASTQ files for accession: $ACCESSION"
#fastq-dump --split-files --gzip --accession "$ACCESSION" --outdir "$RAW_DATA_DIR"


# Step 2: Preprocess the FASTQ files with FastP (Trimming)
for ReadPair in "${RAW_DATA_DIR}/${ACCESSION}"_1.fastq.gz "${RAW_DATA_DIR}/${ACCESSION}"_2.fastq.gz
do
  BASENAME=$(basename "$ReadPair" "_1.fastq.gz")
  echo "Running FastP on $BASENAME"
  fastp \
    -i "${RAW_DATA_DIR}/${ACCESSION}_1.fastq.gz" \
    -o "${FILTERED_DATA_DIR}/Filtered_${ACCESSION}_1.fastq.gz" \
    -I "${RAW_DATA_DIR}/${ACCESSION}_2.fastq.gz" \
    -O "${FILTERED_DATA_DIR}/Filtered_${ACCESSION}_2.fastq.gz" \
    --trim_front1 10 \
    --trim_front2 10

  # Verify FastP ran successfully before moving reports
  if [ -f "fastp.html" ] && [ -f "fastp.json" ]; then
    # Save QC reports for each pair
    mv fastp.html "$REPORT_DIR/${ACCESSION}.html"
    mv fastp.json "$REPORT_DIR/${ACCESSION}.json"
  else
    echo "FastP did not produce expected output files for ${ACCESSION}. Skipping report move."
  fi
done


# Step 3: Check for existing index files with .fasta extension
INDEX_PREFIX="${HOST_REF_PATH}"  # Use the full path including .fasta

if [ -n "$INDEX_DIR" ] && [ -f "${INDEX_DIR}/blackcap_ref_reg.fasta.bwt" ]; then
    echo "Using existing BWA index in $INDEX_DIR"
    INDEX_PREFIX="${INDEX_DIR}/blackcap_ref_reg.fasta"
elif [ -f "${INDEX_PREFIX}.bwt" ] && [ -f "${INDEX_PREFIX}.pac" ] && [ -f "${INDEX_PREFIX}.ann" ] && [ -f "${INDEX_PREFIX}.amb" ] && [ -f "${INDEX_PREFIX}.sa" ]; then
    echo "Using existing BWA index alongside the reference genome"
else
    echo "Indexing reference genome with BWA"
    bwa index "$HOST_REF_PATH"
    if [ $? -ne 0 ]; then
        echo "Error: BWA index creation failed."
        exit 1
    fi
fi

# Step 4: Align filtered reads to the reference genome using BWA
echo "Aligning reads for ${ACCESSION} to the reference genome"
bwa mem -t ${SLURM_CPUS_PER_TASK} "$INDEX_PREFIX" \
    "${FILTERED_DATA_DIR}/Filtered_${ACCESSION}_1.fastq.gz" \
    "${FILTERED_DATA_DIR}/Filtered_${ACCESSION}_2.fastq.gz" | \
    samtools view -bS - > "${RAW_DATA_DIR}/${ACCESSION}.bam"

# Verify if alignment succeeded
if [ $? -ne 0 ]; then
    echo "Error: Alignment with BWA failed for ${ACCESSION}."
    exit 1
fi

# Step 5: Sort and index the final BAM file
echo "Sorting and indexing BAM file for ${ACCESSION}"
samtools sort "${RAW_DATA_DIR}/${ACCESSION}.bam" -o "${RAW_DATA_DIR}/${ACCESSION}.sorted.bam"
samtools index "${RAW_DATA_DIR}/${ACCESSION}.sorted.bam"

# Verify if sorting and indexing succeeded
if [ $? -ne 0 ]; then
    echo "Error: Sorting or indexing the BAM file failed for ${ACCESSION}."
    exit 1
fi

# Optional: Cleanup raw and filtered FASTQ files to save space
# rm ${RAW_DATA_DIR}/${ACCESSION}_*.fastq.gz
# rm ${FILTERED_DATA_DIR}/Filtered_${ACCESSION}_*.fastq.gz

echo "Processing for accession ${ACCESSION} completed"
