#!/bin/bash
# Script to download and process NCBI genomes for various taxonomic groups

# Ensure the script exits if there is an error or if any variables are unset
set -eu

# Input arguments
PARAMETERS_FILE=$1
OUTPUT_PATHS=$2

# Extract the download directory and conda environment path from the config file
DOWNLOAD_DIR=$(grep 'download_dir' "$PARAMETERS_FILE" | cut -d'=' -f2)
NCBI_DATASETS_CONDA=$(grep 'ncbi_datasets_conda' "$PARAMETERS_FILE" | cut -d'=' -f2)

# Ensure that the download directory is set
if [ -z "$DOWNLOAD_DIR" ]; then
    echo "Error: download_dir not defined in config."
    exit 1
fi

# Function to log final clean file paths to the output file
log_file() {
    echo "$1" >> "$OUTPUT_PATHS"
}

# NCBI directory within the specified download directory
NCBI_DIR="$DOWNLOAD_DIR/ncbi"
mkdir -p "$NCBI_DIR"

# Load the required conda environment
module load Anaconda3
source activate "$NCBI_DATASETS_CONDA"

# Taxa list to download
declare -A TAXA=(
    ["apicomplexa"]="wgs_apicomplexa.fasta"
    ["cercozoa"]="wgs_cercozoa.fasta"
    ["euglenozoa"]="wgs_euglenozoa.fasta"
    ["fornicata"]="wgs_fornicata.fasta"
    ["parabasalia"]="wgs_parabasalia.fasta"
)

# Loop over each taxon and download genome data
for taxon in "${!TAXA[@]}"; do
    echo "Downloading genome data for taxon: $taxon"

    # Define directory for the current taxon
    TAXON_DIR="$NCBI_DIR/$taxon"
    mkdir -p "$TAXON_DIR"

    # Download the genome data for the current taxon using the datasets tool
    datasets download genome taxon "$taxon" --dehydrated --filename "$TAXON_DIR/$taxon.zip"

    # Unzip the downloaded file
    unzip "$TAXON_DIR/$taxon.zip" -d "$TAXON_DIR"

    # Move and combine all .fna files into a single file
    find "$TAXON_DIR" -name "*.fna" -exec cat {} + > "$NCBI_DIR/${TAXA[$taxon]}"

    # Log the final FASTA file
    log_file "$NCBI_DIR/${TAXA[$taxon]}"

    # Clean up the unzipped directory and intermediate files
    rm -r "$TAXON_DIR"/*
done

# Deactivate the conda environment
conda deactivate

echo "NCBI genome downloads and processing complete. Final files logged in $OUTPUT_PATHS."
