#!/bin/bash

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --time=3:00:00
#SBATCH --mem 200G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL

# Script to classify sequences using Kraken2.
# Usage:
#   sbatch --array=1-N%M kraken_classify.sh <parameters_file>

set -eu

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <parameters_file>"
    exit 1
fi

PARAMETERS_FILE=$1

# Load parameters from the file
declare -A PARAMS
while IFS="=" read -r key value; do
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$PARAMETERS_FILE"

# Retrieve the required parameters
KRAKEN_CUSTOM=${PARAMS["kraken_custom"]}
OUTPUT_DIR=${PARAMS["output_dir"]}

if [ -z "$KRAKEN_CUSTOM" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Load necessary modules
module purge
module load Kraken2/2.1.1-gompi-2020b

# Directory for Kraken2 output
KRAKEN_OUTPUT="$OUTPUT_DIR/kraken_output"
KRAKEN_REPORT="$OUTPUT_DIR/kraken_reports"

mkdir -p "$KRAKEN_OUTPUT"
mkdir -p "$KRAKEN_REPORT"

# List all assembled contig files (final.contigs.fa)
CONTIG_FILES=("$OUTPUT_DIR/combined_sequences/"*.sorted.combined_sequences.fasta)
NUM_CONTIGS=${#CONTIG_FILES[@]}

if [ "$NUM_CONTIGS" -eq 0 ]; then
    echo "Error: No contig files found."
    exit 1
fi

# Use SLURM_ARRAY_TASK_ID to pick the specific contig file
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
contig_file=${CONTIG_FILES[$INDEX]}

if [ -z "$contig_file" ]; then
    echo "Error: No contig file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Extract sample name from contig file path
sample=$(basename "$contig_file" .sorted.combined_sequences.fasta)

# Run Kraken2 classification
kraken2 --db "$KRAKEN_CUSTOM" \
    --threads 16 \
    --use-names \
    --output "$KRAKEN_OUTPUT/$sample.kraken_output" \
    --report "$KRAKEN_REPORT/$sample.kraken_report" \
    "$contig_file"
