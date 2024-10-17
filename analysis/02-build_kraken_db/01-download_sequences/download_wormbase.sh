#!/bin/bash
# Script to download and process WormBase ParaSite genomes

# Ensure the script exits if there is an error or if any variables are unset
set -eu

# Check if the correct number of arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <parameters_file> <output_paths>"
    exit 1
fi

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

# WormBase directory within the specified download directory
WORMBASE_DIR="$DOWNLOAD_DIR/wormbase"
mkdir -p "$WORMBASE_DIR"

# Step 1: Download the WormBase ParaSite genome data (most recent release)
echo "Downloading WormBase genomes..."
wget -P "$WORMBASE_DIR" -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/ > "$WORMBASE_DIR/wget_log.txt" 2>&1

# Check if the download was successful
echo "Checking download log for any errors..."
cat "$WORMBASE_DIR/wget_log.txt"

# Step 2: Decompress the downloaded genome files
echo "Decompressing WormBase genome files..."

# Use find to locate all .fa.gz files in subdirectories and decompress them
find "$WORMBASE_DIR" -name "*.fa.gz" -exec gunzip -k {} \;

# Step 3: Process headers and concatenate files
echo "Processing headers and concatenating files..."
find "$WORMBASE_DIR" -name "*.fa" | while read decompressed_file; do
    filename=$(basename "$decompressed_file" .fa)
    bioproject=$(echo "$filename" | grep -oP "PRJ[A-Z0-9]+")
    sed -i "s/^>/>$bioproject|/" "$decompressed_file"
done

# Concatenate all .fa files into one final FASTA file
FINAL_FILE="$WORMBASE_DIR/wormbase.fasta"
find "$WORMBASE_DIR" -name "*.fa" -exec cat {} + > "$FINAL_FILE"

# Step 4: Log the final cleaned FASTA file to the output paths file
echo "$FINAL_FILE" >> "$OUTPUT_PATHS"

echo "WormBase download and processing complete. Final file saved at $FINAL_FILE."
