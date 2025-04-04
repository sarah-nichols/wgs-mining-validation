#!/bin/bash
#SBATCH --job-name=download_all
#SBATCH --output=logs/all_%j.out
#SBATCH --error=logs/all_%j.err
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --partition=medium

# Load the config
PARAMETERS_FILE=$1
DOWNLOAD_DIR=$(grep 'download_dir' "$PARAMETERS_FILE" | cut -d'=' -f2)

# Initialize the output paths file (must be in the download directory)
OUTPUT_PATHS="$DOWNLOAD_DIR/downloaded_files.txt"
> "$OUTPUT_PATHS"  # Create or clear the file before starting

# Ensure that the download directory is set
if [ -z "$DOWNLOAD_DIR" ]; then
    echo "Error: download_dir not defined in config."
    exit 1
fi


# Run each download script, passing the parameters file and the output paths file
sbatch /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/analysis/02-build_kraken_db/01-download_sequences/download_eupath.sh "$PARAMETERS_FILE" "$OUTPUT_PATHS"
sbatch /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/analysis/02-build_kraken_db/01-download_sequences/download_ncbi_genomes.sh "$PARAMETERS_FILE" "$OUTPUT_PATHS"
sbatch /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/analysis/02-build_kraken_db/01-download_sequences/download_pr2.sh "$PARAMETERS_FILE" "$OUTPUT_PATHS"
#sbatch /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/analysis/02-build_kraken_db/01-download_sequences/download_ncbi_entrez.sh "$PARAMETERS_FILE" "$OUTPUT_PATHS"
#sbatch /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/analysis/02-build_kraken_db/01-download_sequences/download_wormbase.sh "$PARAMETERS_FILE" "$OUTPUT_PATHS"

