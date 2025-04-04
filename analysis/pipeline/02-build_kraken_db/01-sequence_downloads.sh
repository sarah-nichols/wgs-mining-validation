#!/bin/bash

#SBATCH --job-name=download_dbs
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --mem=80G
#SBATCH --mail-type=ALL

set -eu  # Exit on error, and treat unset variables as an error
set -x   # Enable debugging mode to trace commands

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <parameters_file> [--mail-user=email@example.com]"
    exit 1
fi

# Assign input arguments to variables
PARAMETERS_FILE=$1

# Load parameters from the parameters file
declare -A PARAMS
while IFS="=" read -r key value; do
    # Only assign values if both key and value are not empty
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$PARAMETERS_FILE"

# Retrieve the required parameters from the file
DOWNLOAD_DIR=${PARAMS["download_dir"]}
DATABASE_LIST=${PARAMS["database_list"]}
NCBI_DATASETS_CONDA=${PARAMS["ncbi_datasets_conda"]}
TAXONKIT_CONDA=${PARAMS["taxonkit_conda"]}
ENTREZ_CONDA=${PARAMS["entrez_conda"]}
ENTREZ_PY_SCRIPT=${PARAMS["entrez_py_script"]}

# Validate that all required parameters are provided
if [ -z "$DOWNLOAD_DIR" ] || [ -z "$DATABASE_LIST" ] || [ -z "$NCBI_DATASETS_CONDA" ] || [ -z "$TAXONKIT_CONDA" ] || [ -z "$ENTREZ_CONDA" ] || [ -z "$ENTREZ_PY_SCRIPT" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Define the output paths file
OUTPUT_PATHS="$DOWNLOAD_DIR/downloaded_files.txt"

# Create necessary download directories and output file for paths
mkdir -p "$DOWNLOAD_DIR"
> "$OUTPUT_PATHS"  # Clear the file if it exists

# Function to log file paths
log_file() {
    echo "$1" >> "$OUTPUT_PATHS"
}


# Function to download and process EuPathDB
download_eupath() {
    local EUPATH_DIR="$DOWNLOAD_DIR/eupath"
    mkdir -p "$EUPATH_DIR"

    # Download EuPathDB datasets
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/AmoebaDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/CryptoDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/GiardiaDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/MicrosporidiaDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PiroplasmaDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PlasmoDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/ToxoDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TrichDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TriTrypDB46.tgz
    wget -P "$EUPATH_DIR" ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/seqid2taxid.map

    # Decompress and process files
    for file in "$EUPATH_DIR"/*.tgz; do
        tar -xvf "$file" -C "$EUPATH_DIR"
    done

    # Combine all .fna files into a single FASTA file
    shopt -s globstar
    cat "$EUPATH_DIR"/**/*.fna > "$EUPATH_DIR/eupath.fasta"
    log_file "$EUPATH_DIR/eupath.fasta"

    # Clean up by removing decompressed files except the final combined FASTA
    find "$EUPATH_DIR" -type f ! -name 'eupath.fasta' -exec rm -f {} +

    # Format headers in the FASTA file
    awk '/^>/ {if(seq && ok) {print header ORS seq} header=$0; seq=""; ok = /^>[A-Z0-9]+(\.[0-9]+)? \|/} /^[^>]/ {seq = seq ? seq ORS $0 : $0; ok = ok && !/^N+$/ && !/^x+$/} END {if(seq && ok) print header ORS seq}' "$EUPATH_DIR/eupath.fasta" > "$EUPATH_DIR/eupath_headers.fasta"
    log_file "$EUPATH_DIR/eupath_headers.fasta"
}



# Function to download and process NCBI genomes
download_ncbi_genomes() {
    local NCBI_DIR="$DOWNLOAD_DIR/ncbi"
    mkdir -p "$NCBI_DIR"

    module purge
    module load Anaconda3/2024.02-1
    source activate "$NCBI_DATASETS_CONDA"

    # Download and process sequences for each taxon
    declare -A TAXA=(
        ["apicomplexa"]="wgs_apicomplexa.fasta"
        ["cercozoa"]="wgs_cercozoa.fasta"
        ["euglenozoa"]="wgs_euglenozoa.fasta"
        ["fornicata"]="wgs_fornicata.fasta"
        ["parabasalia"]="wgs_parabasalia.fasta"
    )

    for taxon in "${!TAXA[@]}"; do
    echo "Processing $taxon..."

    # Download and unzip
    datasets download genome taxon "$taxon" || true  # Continue if this step fails
    unzip ncbi_dataset.zip -d "$NCBI_DIR/$taxon" || true  # Continue if the unzip step fails

    # Debugging: Check files after decompression
    echo "Files in $NCBI_DIR/$taxon after decompression:"
    ls "$NCBI_DIR/$taxon" || true  # Continue if this step fails

    # Move .fna files
    find "$NCBI_DIR/$taxon" -type f -name "*.fna" -exec mv {} "$NCBI_DIR/$taxon" \; || true  # Continue if no .fna files are found

    # Concatenate the .fna files
    echo "Concatenating .fna files for $taxon"
    cat "$NCBI_DIR/$taxon/*.fna" > "$NCBI_DIR/${TAXA[$taxon]}" || true  # Continue even if concatenation fails
    log_file "$NCBI_DIR/${TAXA[$taxon]}" || true  # Continue if logging fails

    # Clean up subdirectories
    echo "Cleaning up directories for $taxon"
    rm -r "$NCBI_DIR/$taxon/" || true  # Continue if the cleanup fails
done

    conda deactivate
}

set -e


# Function to download and process PR2 database
download_pr2() {
    local PR2_DIR="$DOWNLOAD_DIR/PR2"
    mkdir -p "$PR2_DIR"

    wget -O "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
    gzip -d "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
    log_file "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta"

    # Filter sequences
    awk 'BEGIN{IGNORECASE=1} /^>/{printit=1; for (i=1;i<=NF;i++) {if ($i ~ /16S_rRNA|Fungi|uncultured|Arthropoda|Annelida|Craniata|Streptophyta|Chlorophyta/) {printit=0}} } {if (printit) print}' "$PR2_DIR/pr2_version_5.0.0_SSU_taxo_long.fasta" > "$PR2_DIR/pr2_clean.fasta"
    log_file "$PR2_DIR/pr2_clean.fasta"
}

# Function to download and process WormBase
download_wormbase() {
    local WORMBASE_DIR="$DOWNLOAD_DIR/wormbase"
    mkdir -p "$WORMBASE_DIR"

    wget -P "$WORMBASE_DIR" -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

    shopt -s globstar
    set +e
    for file in "$WORMBASE_DIR/**/*.fa.gz"; do
        gunzip -k "$file"
        decompressed_file="${file%.gz}"
        filename=$(basename "$decompressed_file" .fa)
        bioproject=$(echo ${filename} | grep -oP "PRJ[A-Z]+[0-9]+")
        sed -i "s/^>/>$bioproject|/" "$decompressed_file"
    done
    set -e

    cat "$WORMBASE_DIR/**/*.fa" > "$WORMBASE_DIR/wormbase.fasta"
    log_file "$WORMBASE_DIR/wormbase.fasta"
    rm -r "$WORMBASE_DIR/ftp.ebi.ac.uk"
}

# Step 1: Generate the Taxonomy IDs list
echo "Generating Taxonomy IDs list..."

module load Anaconda3/2024.02-1
source activate "$TAXKIT_ENV"

taxonkit list \
    --ids 5794,207245,5719,33682,136419 \
    --indent "" --data-dir "$TAXDUMP_DIR" \
    --show-name | grep -v -E 'uncultured|environmental' | awk '{print $1}' > "$TAXID_FILE"

conda deactivate

echo "Taxonomy ID list generated at $TAXID_FILE"

# Step 2: Run the Python script using the generated taxids
echo "Running the Python sequence download script..."

source activate "$ENTREZ_ENV"
python "$ENTREZ_PY_SCRIPT" "$PARAMETERS_FILE"
conda deactivate

echo "Sequence download completed."

# Function to run NCBI sequence downloads via Entrez
download_ncbi_entrez() {
    local NCBI_ENTREZ_DIR="$DOWNLOAD_DIR/ncbi_entrez"
    mkdir -p "$NCBI_ENTREZ_DIR"

    module purge
    module load Anaconda3/2024.02-1

    # Activate the conda environment for taxonkit
    source activate "$TAXONKIT_CONDA"

    taxonkit list --ids 5794,207245,5719,33682,136419 --indent "" \
        --data-dir /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/taxid_map/taxdump \
        --show-name | grep -v -E 'uncultured|environmental' |
        awk '{print $1}' > "$NCBI_ENTREZ_DIR/taxids_include.txt"
    log_file "$NCBI_ENTREZ_DIR/taxids_include.txt"

    conda deactivate

    # Activate conda environment for Python Entrez script
    source activate "$ENTREZ_CONDA"
    python "$ENTREZ_PY_SCRIPT"
    conda deactivate
}


# Process the database list from DATABASE_LIST
IFS=',' read -ra DBS <<< "$DATABASE_LIST"
for db in "${DBS[@]}"; do
    case "$db" in
        eupath)
            download_eupath
            ;;
        ncbi_genomes)
            download_ncbi_genomes
            ;;
        pr2)
            download_pr2
            ;;
        wormbase)
            download_wormbase
            ;;
        ncbi_entrez)
            download_ncbi_entrez
            ;;
        *)
            echo "Unknown database: $db"
            ;;
    esac
done

echo "Download and processing complete."
echo "All paths to generated files have been saved to $OUTPUT_PATHS."