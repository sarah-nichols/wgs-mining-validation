#!/bin/bash
# Script to download and process EuPathDB datasets

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

# EuPathDB directory within the specified download directory
EUPATH_DIR="$DOWNLOAD_DIR/eupath"
mkdir -p "$EUPATH_DIR"

# Download EuPathDB datasets
echo "Downloading EuPathDB datasets..."
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/AmoebaDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/CryptoDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/GiardiaDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/MicrosporidiaDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PiroplasmaDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PlasmoDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/ToxoDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TrichDB46.tgz
wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TriTrypDB46.tgz

# Uncompress the downloaded datasets
echo "Decompressing EuPathDB datasets..."
for file in "$EUPATH_DIR"/*.tgz; do
    tar -xvf "$file" -C "$EUPATH_DIR"
done

# Count and log number of headers in each .fna file
echo "Counting headers in individual .fna files..."
for fna_file in "$EUPATH_DIR"/**/*.fna; do
    count=$(grep -c "^>" "$fna_file")
    echo "File: $fna_file contains $count sequences."
done

# Combine all .fna files into a single FASTA file
echo "Combining .fna files..."
cat "$EUPATH_DIR"/**/*.fna > "$EUPATH_DIR/eupath_combined.fasta"

# Count and log number of headers in the combined FASTA file
combined_count=$(grep -c "^>" "$EUPATH_DIR/eupath_combined.fasta")
echo "Combined FASTA file contains $combined_count sequences."

# Prepare file for removed headers
REMOVED_HEADERS_FILE="$EUPATH_DIR/removed_headers.txt"
> "$REMOVED_HEADERS_FILE" # Create an empty file to store removed headers

# Clean headers in the combined FASTA file (modify as needed for your specific header format)
echo "Cleaning FASTA headers..."
awk '
    /^>/ { 
        if(seq && ok) {
            print header ORS seq
        } 
        header=$0; seq=""; ok = /^>[A-Z0-9]+(\.[0-9]+)? \|/ 
        if(!ok) {
            print header >> "'$REMOVED_HEADERS_FILE'"
        }
    } 
    /^[^>]/ { 
        seq = seq ? seq ORS $0 : $0; 
        ok = ok && !/^N+$/ && !/^x+$/ 
    } 
    END { 
        if(seq && ok) {
            print header ORS seq
        }
    }' "$EUPATH_DIR/eupath_combined.fasta" > "$EUPATH_DIR/eupath_clean_headers.fasta"

# Count and log number of headers in the cleaned FASTA file
cleaned_count=$(grep -c "^>" "$EUPATH_DIR/eupath_clean_headers.fasta")
echo "Cleaned FASTA file contains $cleaned_count sequences."

# Count and log number of removed headers
removed_count=$(grep -c "^>" "$REMOVED_HEADERS_FILE")
echo "Removed $removed_count headers. They are stored in $REMOVED_HEADERS_FILE."

# Log the final cleaned FASTA file
log_file "$EUPATH_DIR/eupath_clean_headers.fasta"

# Remove intermediate files to save space, keeping only the final clean FASTA file and cleaned headers file
echo "Cleaning up intermediate files..."
find "$EUPATH_DIR" -type f ! -name 'eupath_clean_headers.fasta' ! -name 'eupath_clean_headers.txt' -exec rm -f {} +

echo "EuPathDB download and processing complete."
echo "Final file saved to $EUPATH_DIR/eupath_clean_headers.fasta"