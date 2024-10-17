#!/bin/sh

#SBATCH --job-name=remove_irreg_exp
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL

# Script to remove commas from an input file and save the result to an output file
# Usage:
#   sbatch this_script.sh <input_file> <output_file> [email]
# 
# Arguments:
#   <input_file>  - Path to the input FASTA file where commas will be removed
#   <output_file> - Path to the output file where the processed data will be saved
#   [email]       - (Optional) Email address to receive job notifications

set -eu

# Check if the correct number of arguments are provided
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <input_file> <output_file> [email]"
    exit 1
fi

# Assign input arguments to variables
INPUT_FILE=$1
OUTPUT_FILE=$2
EMAIL=${3:-}

# Set the email if provided
if [ -n "$EMAIL" ]; then
    # Add the --mail-user directive dynamically if an email is provided
    sbatch --mail-user="$EMAIL"
fi

# Display the paths of the input and output files
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"

# Run the sed command to remove commas from the input file and save the result in the output file
sed 's/,//' "$INPUT_FILE" > "$OUTPUT_FILE"
