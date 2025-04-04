#!/bin/bash

#SBATCH --job-name=build_host_ref
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 80G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL

# Script to build host reference using GATK's BwaMemIndexImageCreator and PathSeqBuildKmers tools.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file> - Path to the text file containing paths to the input reference file and output files
#   [email]           - (Optional) Email address to receive job notifications

# Load necessary modules
module purge
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
INPUT_REG=${PARAMS["input_reg"]}
OUTPUT_IMG=${PARAMS["output_img"]}
OUTPUT_HSS=${PARAMS["output_hss"]}

# Validate that all required parameters are provided
if [ -z "$INPUT_REG" ] || [ -z "$OUTPUT_IMG" ] || [ -z "$OUTPUT_HSS" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

echo "Building host reference..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
    -I  "$INPUT_REG" \
    -O  "$OUTPUT_IMG"

echo "Output file: $OUTPUT_IMG"

echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$INPUT_REG" \
   --output "$OUTPUT_HSS" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "PathSeqBuildKmers completed"
