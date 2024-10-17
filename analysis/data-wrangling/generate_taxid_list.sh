#!/bin/bash

# Set the paths for the taxonkit Conda environment, taxdump directory, and output file
export TAXONKIT_CONDA="/data/biol-bird-parasites/sann7416/conda_environments/taxonkit_env"
export TAXDUMP_DIR="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/taxid_map/taxdump"
export TAXID_FILE="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/taxids_not_interested.txt"

# Load Anaconda3 and activate the taxonkit environment
module load Anaconda3
source activate "$TAXONKIT_CONDA"

# Check if taxdump directory exists
if [ ! -d "$TAXDUMP_DIR" ]; then
    echo "Error: TAXDUMP_DIR does not exist."
    exit 1
fi

# List taxonomy using taxonkit and filter for Bacteria (TaxID: 2) and Vertebrata (TaxID: 7742)
taxonkit list --ids 2,7742 --indent "" --data-dir "$TAXDUMP_DIR" --show-name \
    | grep -v -E 'uncultured|environmental' \
    | awk '{print $1}' > "$TAXID_FILE"

# Check if the output file was generated and contains data
if [ -s "$TAXID_FILE" ]; then
    echo "Bacteria (TaxID 2) and Vertebrate (TaxID 7742) Taxonomy IDs list saved to $TAXID_FILE"
else
    echo "Warning: No data written to $TAXID_FILE."
fi

# Deactivate the conda environment
conda deactivate
