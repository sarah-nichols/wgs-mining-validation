#!/bin/sh

#SBATCH --job-name=random_sub_sample
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --time=5:00:00
#SBATCH --mem=30G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --array=1-11  # Replace N with the number of BAM files
#SBATCH --mail-type=ALL

# Load SAMtools
module purge
module load SAMtools/1.14-GCC-11.2.0

# Check if the correct number of arguments are provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <config_file.txt> [--mail-user=email@example.com]"
    exit 1
fi

# Assign input arguments to variables
CONFIG_FILE=$1
EMAIL=${2:-}

# Set the email if provided
if [ -n "$EMAIL" ]; then
    sbatch --mail-user="$EMAIL"
fi

# Load parameters from the configuration file
declare -A PARAMS
while IFS="=" read -r key value; do
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$CONFIG_FILE"

# Retrieve parameters from the configuration file
BAM_DIR=${PARAMS["bam_dir"]}
OUTPUT_DIR=${PARAMS["output_dir"]}
COVERAGES=${PARAMS["coverages"]}

# Validate that all required parameters are provided
if [ -z "$BAM_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$COVERAGES" ]; then
    echo "Error: Missing required parameters in the configuration file."
    exit 1
fi

# Create the output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Convert coverages string to an array
IFS=',' read -ra COVERAGE_ARRAY <<< "$COVERAGES"

# Get the BAM file corresponding to the current SLURM_ARRAY_TASK_ID
BAM_FILES=("$BAM_DIR"/*.sorted.bam)
BAM_FILE=${BAM_FILES[$((SLURM_ARRAY_TASK_ID - 1))]}

# Ensure the BAM file exists
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file does not exist for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Get the base name of the BAM file
bam_base=$(basename "$BAM_FILE" .sorted.bam)
echo "Processing $bam_base"

# Set genome size in base pairs
GENOME_SIZE=1090000000  # 1.09 Gbp

# Define an average read length (assuming 150 bp)
average_read_length=150

# Get the total number of reads in the BAM file
total_reads=$(samtools idxstats "$BAM_FILE" | awk '{s+=$3} END {print s}')
echo "Total reads in $bam_base: $total_reads"

# Calculate the total number of bases covered by all reads
total_bases=$(echo "$total_reads * $average_read_length" | bc)  
current_coverage=$(echo "$total_bases / $GENOME_SIZE" | bc -l)
echo "Current coverage: ${current_coverage}X"

# Calculate maximum coverage
max_coverage=$(echo "$total_bases / $GENOME_SIZE" | bc -l)
echo "Maximum coverage possible for $bam_base: ${max_coverage}X"

# Loop over each desired coverage
for coverage in "${COVERAGE_ARRAY[@]}"; do
    # Skip if the desired coverage is greater than the max coverage
    if [ "$(echo "$coverage > $max_coverage" | bc)" -eq 1 ]; then
        echo "Skipping $coverage X as it exceeds the maximum coverage of ${max_coverage}X"
        continue
    fi

    # Calculate the number of reads required for the desired coverage
    target_reads=$(echo "$coverage * $GENOME_SIZE / $average_read_length" | bc)

    # Calculate the subsample fraction
    subsample_fraction=$(echo "$target_reads / $total_reads" | bc -l)

    # Ensure the subsample fraction is valid
    if [ -z "$subsample_fraction" ] || (( $(echo "$subsample_fraction <= 0" | bc -l) )); then
        echo "Skipping subsampling as the subsample fraction is invalid: $subsample_fraction"
        continue
    fi

    echo "Subsampling $bam_base to approximately ${coverage}X coverage"

    # Create a subsampled BAM file
    output_bam="${OUTPUT_DIR}/${bam_base}_X${coverage}.bam"
    samtools view -s "$subsample_fraction" "$BAM_FILE" -o "$output_bam"
done
