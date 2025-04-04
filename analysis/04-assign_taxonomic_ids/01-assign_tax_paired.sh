#!/bin/bash

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --time=3:00:00
#SBATCH --mem 180G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL

# Script to process BAM files, convert them to FASTQ, and classify sequences using Kraken2.
# Usage:
#   sbatch --array=1-N%M kraken_assign.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file>    - Path to the text file containing paths to BAM files, Kraken2 database, and output directory.
#   [email]              - (Optional) Email address to receive job notifications.
#
# The SLURM array length (N) and concurrency limit (M) should be specified by the user at job submission.

set -eu

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <parameters_file> [--mail-user=email@example.com]"
    exit 1
fi

# Assign input arguments to variables
PARAMETERS_FILE=$1

# Optionally accept email as a command-line argument (if provided, it will be used in the job submission)
EMAIL=${2:-}

# Load parameters from the file
declare -A PARAMS
while IFS="=" read -r key value; do
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$PARAMETERS_FILE"

# Retrieve the required parameters from the file
BAM_PATH=${PARAMS["bam_path"]}
KRAKEN_CUSTOM=${PARAMS["kraken_custom"]}
OUTPUT_DIR=${PARAMS["output_dir"]}

# Validate that all required parameters are provided
if [ -z "$BAM_PATH" ] || [ -z "$KRAKEN_CUSTOM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Load necessary modules
module purge
module load Kraken2/2.1.1-gompi-2020b
module load SAMtools/1.16.1-GCC-11.3.0
module load SeqKit/2.2.0

# Create necessary output directories
PAIRED_F_FASTQS="$OUTPUT_DIR/host_filtered/paired_fastqs_forward"
PAIRED_R_FASTQS="$OUTPUT_DIR/host_filtered/paired_fastqs_reverse"
KRAKEN_OUTPUT_PAIRED="$OUTPUT_DIR/kraken_output/output"
KRAKEN_REPORT_PAIRED="$OUTPUT_DIR/kraken_output/reports"

mkdir -p "$PAIRED_F_FASTQS"
mkdir -p "$PAIRED_R_FASTQS"
mkdir -p "$KRAKEN_OUTPUT_PAIRED"
mkdir -p "$KRAKEN_REPORT_PAIRED"

# List all BAM files in the directory
BAM_FILES=("$BAM_PATH"/paired_*.bam)
NUM_FILES=${#BAM_FILES[@]}

# Validate that there are BAM files to process
if [ "$NUM_FILES" -eq 0 ]; then
    echo "Error: No BAM files found in the directory."
    exit 1
fi

# Use SLURM_ARRAY_TASK_ID to pick a specific file
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
bam_file=${BAM_FILES[$INDEX]}

# Ensure bam_file is not empty or unset
if [ -z "$bam_file" ]; then
    echo "Error: No BAM file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Construct output FASTQ file names
forward_fastq_file="$PAIRED_F_FASTQS/$(basename "${bam_file%.bam}").1.fastq"
reverse_fastq_file="$PAIRED_R_FASTQS/$(basename "${bam_file%.bam}").2.fastq"

# Convert BAM to FASTQ using samtools
samtools fastq -1 "$forward_fastq_file" -2 "$reverse_fastq_file" -n "$bam_file"

sample=$(basename "${forward_fastq_file%.1.fastq}.merged" .fastq)

# Run Kraken2 for classification
kraken2 --db "$KRAKEN_CUSTOM" \
    --threads 8 \
    --paired \
    --use-names \
    --output "$KRAKEN_OUTPUT_PAIRED/$sample.stats.report" \
    --report "$KRAKEN_REPORT_PAIRED/$sample.report" \
    "${forward_fastq_file}" "${reverse_fastq_file}"
