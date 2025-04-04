#!/bin/bash
# Script to download and process NCBI sequences using Entrez, with taxon filtering via TaxonKit

# Ensure the script exits if there is an error or if any variables are unset
set -eu

# Input arguments
PARAMETERS_FILE=$1
OUTPUT_PATHS=$2

# Extract necessary parameters from the config file
DOWNLOAD_DIR=$(grep 'download_dir' "$PARAMETERS_FILE" | cut -d'=' -f2)
ENTREZ_CONDA=$(grep 'entrez_conda' "$PARAMETERS_FILE" | cut -d'=' -f2)
TAXONKIT_CONDA=$(grep 'taxonkit_conda' "$PARAMETERS_FILE" | cut -d'=' -f2)
TAXDUMP_DIR=$(grep 'taxdump_dir' "$PARAMETERS_FILE" | cut -d'=' -f2)
TAXID_FILE=$(grep 'taxid_file' "$PARAMETERS_FILE" | cut -d'=' -f2)
ENTREZ_PY_SCRIPT=$(grep 'entrez_py_script' "$PARAMETERS_FILE" | cut -d'=' -f2)

# Ensure that the download directory is set
if [ -z "$DOWNLOAD_DIR" ]; then
    echo "Error: download_dir not defined in config."
    exit 1
fi

# Function to log final clean file paths to the output file
log_file() {
    echo "$1" >> "$OUTPUT_PATHS"
}

# NCBI Entrez directory within the specified download directory
NCBI_ENTREZ_DIR="$DOWNLOAD_DIR/ncbi_entrez"
mkdir -p "$NCBI_ENTREZ_DIR"

# Step 1: Generate the list of Taxonomy IDs using TaxonKit
echo "Generating Taxonomy IDs list with TaxonKit..."

module load Anaconda3
source activate "$TAXONKIT_CONDA"

# Define the directory for taxdump files
TAXDUMP_DIR="$NCBI_ENTREZ_DIR/taxid_map/taxdump"
mkdir -p "$TAXDUMP_DIR"

# Download and extract taxdump.tar.gz
wget -P "$TAXDUMP_DIR" ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf "$TAXDUMP_DIR/taxdump.tar.gz" -C "$TAXDUMP_DIR"
rm "$TAXDUMP_DIR/taxdump.tar.gz"

# Ensure the required files are in the correct directory
if [ ! -f "$TAXDUMP_DIR/names.dmp" ] || [ ! -f "$TAXDUMP_DIR/nodes.dmp" ] || [ ! -f "$TAXDUMP_DIR/delnodes.dmp" ] || [ ! -f "$TAXDUMP_DIR/merged.dmp" ]; then
    echo "Error: Required taxonomy data files are missing in $TAXDUMP_DIR"
    exit 1
fi

taxonkit list --ids 5794,207245,5719,33682,136419 \
    --indent "" --data-dir "$TAXDUMP_DIR" --show-name \
    | grep -v -E 'uncultured|environmental' | awk '{print $1}' > "$TAXID_FILE"

echo "Taxonomy IDs list saved to $TAXID_FILE"

conda deactivate

# Step 2: Run the Entrez download script using the generated taxids
echo "Running Entrez download script..."
source activate "$ENTREZ_CONDA"

python "$ENTREZ_PY_SCRIPT" --taxid-file "$TAXID_FILE" --output "$NCBI_ENTREZ_DIR/ncbi_sequences.fasta"

# Log the final FASTA file
log_file "$NCBI_ENTREZ_DIR/ncbi_sequences.fasta"

# Deactivate the Conda environment
conda deactivate

echo "NCBI Entrez downloads complete. Final file logged in $OUTPUT_PATHS."
