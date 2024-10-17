#!/bin/bash

#SBATCH --job-name=download_ref
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 10G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL

# Script to download a reference genome, process it, and output a modified version.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file> - Path to the text file containing paths to the necessary directories and environment variables
#   [email]           - (Optional) Email address to receive job notifications

module purge
module load Anaconda3/2024.02-1

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
NCBI_DATASETS_CONDA=${PARAMS["ncbi_datasets_conda"]}
HOST_REF_DIR=${PARAMS["host_ref_dir"]}

# Validate that all required parameters are provided
if [ -z "$NCBI_DATASETS_CONDA" ] || [ -z "$HOST_REF_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Activate the Conda environment for NCBI datasets
source activate "$NCBI_DATASETS_CONDA"

# Download the genome using the accession number and unzip it into the provided directory
datasets download genome accession GCA_009819655.1
unzip ncbi_dataset.zip -d "$HOST_REF_DIR/blackcap_GCA_009819655.1"

# Deactivate the Conda environment
conda deactivate

# Define the input and output file paths
INPUT_FILE="$HOST_REF_DIR/blackcap_GCA_009819655.1/ncbi_dataset/data/GCA_009819655.1/GCA_009819655.1_bSylAtr1.pri_genomic.fna"
OUTPUT_FILE="$HOST_REF_DIR/blackcap_ref_reg.fasta"

# Echo the input and output files
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"

# Modify the reference genome by removing commas and save it to the output file
sed 's/,//' "$INPUT_FILE" > "$OUTPUT_FILE"

# Remove the downloaded directory after processing
rm -r "$HOST_REF_DIR/blackcap_GCA_009819655.1"
