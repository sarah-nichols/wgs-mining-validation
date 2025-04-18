#!/bin/bash

#SBATCH --job-name=build_kraken_db
#SBATCH --partition long
#SBATCH --time=80:00:00
#SBATCH --mem=200G
#SBATCH --clusters=all
#SBATCH --cpus-per-task=48  # Request 48 CPUs for this job
#SBATCH --mail-type=ALL

# Script to build a custom Kraken2 database.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file>    - Path to the text file containing paths to sequence files and the Kraken database directory
#   [email]              - (Optional) Email address to receive job notifications

set -eu

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
SEQS_TO_ADD=${PARAMS["seqs_to_add"]}
KRAKEN_CUSTOM=${PARAMS["kraken_custom"]}
EUPATH_MAP=${PARAMS["eupath_map"]}

# Validate that all required parameters are provided
if [ -z "$SEQS_TO_ADD" ] || [ -z "$KRAKEN_CUSTOM" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Load Kraken2 module
module purge
module load Kraken2/2.1.2-gompi-2021b

# Download taxonomy and add eupath taxid map
#kraken2-build --download-taxonomy --db "$KRAKEN_CUSTOM" --threads 48
#awk -F'|' '{print $3 "\t" $NF}' "$EUPATH_MAP/seqid2taxid.map" > "$EUPATH_MAP/eupath_accession2taxid.txt"
#cat "$EUPATH_MAP/eupath_accession2taxid.txt" >> "$KRAKEN_CUSTOM/taxonomy/nucl_gb_accession2taxid"

#build the standard and protozoa libraries
#kraken2-build --download-library bacteria --db "$KRAKEN_CUSTOM" --threads 48
#kraken2-build --download-library protozoa --db "$KRAKEN_CUSTOM" --threads 48
#kraken2-build --download-library UniVec --db "$KRAKEN_CUSTOM" --threads 48
#kraken2-build --download-library human --db "$KRAKEN_CUSTOM" --threads 48

# Add sequence files to the Kraken2 library
#while IFS= read -r line; do
#    echo "Adding file: $line to library"
#    kraken2-build --add-to-library "$line" --db "$KRAKEN_CUSTOM" --threads 48 || { echo "kraken2-build failed for $line"; exit 1; }
#done < "$SEQS_TO_ADD"

# Build the Kraken2 database
kraken2-build --build --threads 48 --db "$KRAKEN_CUSTOM" --fast-build

# Clean up the Kraken2 database build files
kraken2-build --clean --threads 48 --db "$KRAKEN_CUSTOM"



