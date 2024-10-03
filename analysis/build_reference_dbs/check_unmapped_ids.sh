#!/bin/bash
#SBATCH --job-name=find_unmapped
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition short
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

:<<'COMMENT'

# File containing the list of accession numbers
accession_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/kraken_custom/unmapped.txt"

# Define the paths to the FASTA files
fasta_files=(
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/metabarcoding_clean.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/wormbase/wormbase.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/wgs_apicomplexa.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/wgs_cercozoa.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/wgs_euglenozoa.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/wgs_fornicata.fasta"
    "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/wgs_parabasalia.fasta"
)

# Output file for results
output_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/unmapped_accession_search_results.txt"
> "$output_file"  # Clear the output file

# Step 1: Read accession numbers from unmapped.txt into an array
mapfile -t accession_list < "$accession_file"

# Function to search for accession in a single FASTA file
search_accessions_in_fasta() {
    local fasta_file="$1"
    echo "Searching in $fasta_file..."
    
    # Loop through each accession and search in the current FASTA file
    for accession in "${accession_list[@]}"; do
        if grep -q "$accession" "$fasta_file"; then
            echo "Found $accession in $fasta_file" >> "$output_file"
        fi
    done
}

# Step 2: Loop through all FASTA files and search
for fasta_file in "${fasta_files[@]}"; do
    if [[ -f "$fasta_file" ]]; then
        search_accessions_in_fasta "$fasta_file"
    fi
done

echo "Search completed. Results saved to $output_file."

COMMENT

# Define the file paths
unmapped_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/kraken_custom/unmapped.txt"  # File containing the list of unmapped accession numbers
seqid2taxid_file="/data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/data/reference_database/eupath_map/seqid2taxid.map"  # Path to seqid2taxid.map
output_file="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/accessions_eupath_map.txt"  # Output file to store results

# Clear the output file
> "$output_file"

# Loop through each accession in the unmapped file and check if it's in seqid2taxid.map
while read -r accession; do
    if grep -q "$accession" "$seqid2taxid_file"; then
        echo "$accession found in seqid2taxid.map" >> "$output_file"
    else
        echo "$accession not found in seqid2taxid.map" >> "$output_file"
    fi
done < "$unmapped_file"

echo "Search completed. Results saved to $output_file."
