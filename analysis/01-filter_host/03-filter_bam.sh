#!/bin/bash

#SBATCH --job-name=filter_host
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --time=130:00:00
#SBATCH --mail-type=ALL

# Script to filter host BAM files using GATK's PathSeqFilterSpark.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file> - Path to the text file containing paths to the sample IDs file, host HSS file, host image file, and output directory
#   [email]           - (Optional) Email address to receive job notifications

module purge
module load SAMtools/1.14-GCC-11.2.0
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

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
SAMPLE_IDS_FILE=${PARAMS["sample_ids_file"]}
HOST_HSS=${PARAMS["host_hss"]}
HOST_IMG=${PARAMS["host_img"]}
OUTPUT_DIR=${PARAMS["output_directory"]}

# Validate that all required parameters are provided
if [ -z "$SAMPLE_IDS_FILE" ] || [ -z "$HOST_HSS" ] || [ -z "$HOST_IMG" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Create necessary output directories
FILTERED_BAMS_DIR="$OUTPUT_DIR/filtered_bams"
FILTER_METRICS_DIR="$OUTPUT_DIR/filter_metrics"

mkdir -p "$FILTERED_BAMS_DIR"
mkdir -p "$FILTER_METRICS_DIR"

# Get the $SLURM_ARRAY_TASK_ID-th line from the sample IDs file
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_IDS_FILE")

BAM_OUT="$FILTERED_BAMS_DIR/${BAM##/*/}"
echo "Copying to $BAM_OUT"
cp "$BAM" "$BAM_OUT"

echo "Editing header with samtools"
BAM_HEADER="$FILTERED_BAMS_DIR/header_${BAM##/*/}"
samtools view -H "$BAM_OUT" | sed 's/,//' | samtools reheader - "$BAM_OUT" > "$BAM_HEADER"
rm "$BAM_OUT"

BAM_FILTERED_PAIRED="$FILTERED_BAMS_DIR/paired_${BAM##/*/}"
echo "Paired output: $BAM_FILTERED_PAIRED"
BAM_FILTERED_UNPAIRED="$FILTERED_BAMS_DIR/unpaired_${BAM##/*/}"
echo "Unpaired output: $BAM_FILTERED_UNPAIRED"

gatk --java-options "-Xmx80G" PathSeqFilterSpark  \
    --input "$BAM_HEADER" \
    --paired-output "$BAM_FILTERED_PAIRED" \
    --unpaired-output "$BAM_FILTERED_UNPAIRED" \
    --min-adapter-length 1 \
    --min-clipped-read-length 60 \
    --kmer-file "$HOST_HSS" \
    --is-host-aligned TRUE \
    --filter-bwa-image "$HOST_IMG" \
    --filter-metrics "$FILTER_METRICS_DIR/metrics_${BAM##/*/}" \
    --bam-partition-size 4000000

rm "$BAM_HEADER"
