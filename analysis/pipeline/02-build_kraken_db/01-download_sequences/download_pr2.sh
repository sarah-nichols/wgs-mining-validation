#!/bin/bash
# Script to download and process PR2 database

# Ensure the script exits if there is an error or if any variables are unset
set -eu

# Input arguments
PARAMETERS_FILE=$1
OUTPUT_PATHS=$2

# Extract the download directory from the config file
DOWNLOAD_DIR=$(grep 'download_dir' "$PARAMETERS_FILE" | cut -d'=' -f2)

# Ensure that the download directory is set
if [ -z "$DOWNLOAD_DIR" ]; then
    echo "Error: download_dir not defined in config."
    exit 1
fi

# Function to log final clean file paths to the output file
log_file() {
    echo "$1" >> "$OUTPUT_PATHS"
}

# PR2 directory within the specified download directory
PR2_DIR="$DOWNLOAD_DIR/PR2"
mkdir -p "$PR2_DIR"

# Step 1: Download the PR2 database (version 5.0.0)
echo "Downloading PR2 database..."
wget -O "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"

# Step 2: Uncompress the downloaded FASTA file
echo "Decompressing PR2 database..."
gzip -d "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"

# Step 3: Filter sequences based on taxonomy and keywords
echo "Filtering PR2 sequences..."
awk 'BEGIN{IGNORECASE=1} 
/^>/{printit=1; 
  for (i=1;i<=NF;i++) {
    if ($i ~ /16S_rRNA|Fungi|uncultured|Arthropoda|Annelida|Craniata|Streptophyta|Chlorophyta/) {
      printit=0
      print $0 > "'"$PR2_DIR/removed_headers.fasta"'"
    }
  }} 
{if (printit) print}' "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta" > "$PR2_DIR/pr2_clean.fasta"

# Step 4: Log the final cleaned FASTA file
log_file "$PR2_DIR/pr2_clean.fasta"

echo "PR2 download and processing complete. Final file saved at $PR2_DIR/pr2_clean.fasta."
echo "Removed headers saved at $PR2_DIR/removed_headers.fasta."