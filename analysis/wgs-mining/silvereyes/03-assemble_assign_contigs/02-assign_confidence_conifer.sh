#!/bin/bash

#SBATCH --job-name=conifer_conf
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 20G
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL

# Script to run Conifer on Kraken2 output files and generate a final report.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file> - Path to the text file containing paths to Kraken2 results directory, Kraken2 custom taxid database, and output directory
#   [email]           - (Optional) Email address to receive job notifications

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
KRAKEN_RESULTS=${PARAMS["kraken_results"]}
KRAKEN_CUSTOM_TAXID=${PARAMS["kraken_custom_taxid"]}
OUTPUT_DIR=${PARAMS["output_directory"]}
CONIFER_SOFTWARE=${PARAMS["conifer_software"]}

# Validate that all required parameters are provided
if [ -z "$KRAKEN_RESULTS" ] || [ -z "$KRAKEN_CUSTOM_TAXID" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$CONIFER_SOFTWARE" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Define the output file for the final Conifer results
CONIFER_OUTPUT="$OUTPUT_DIR/conifer_output.tsv"

# Move to the directory where Conifer is located
cd "$CONIFER_SOFTWARE"

# Create an empty file to store the final output
> "$CONIFER_OUTPUT"

# Print the header line to the output file
echo -e "sample_id\tclassified\tbreakdown\ttaxon_name\tread_length\tbreakdown\tC1\tC2\tCA\tRTL1\tRTL2\tRTLA" > "$CONIFER_OUTPUT"

# Loop through all Kraken2 report files in the directory
for file in "$KRAKEN_RESULTS"/*.sorted.kraken_output; do
  # Run Conifer on the file
  ./conifer --both_scores -i "$file" -d "$KRAKEN_CUSTOM_TAXID" > output.txt

  # Get the filename without the path
  filename=$(basename "$file")

  # Extract the sample ID from the filename
  sample_id=$(echo $filename | awk '{match($0, /^[A-Za-z]{3}[0-9]{1,2}/, arr); print arr[0]}')

  echo "Processing file: $filename"
  echo "Sample ID: $sample_id"

  # Add a column for the sample ID to the output file, skipping the header line
  awk -v id="$sample_id" 'NR>1 {printf("%s\t%s\n", id, $0)}' output.txt > output_with_id.txt

  # Append the output to the final output file
  cat output_with_id.txt >> "$CONIFER_OUTPUT"
done


echo "Conifer processing complete. Results saved to $CONIFER_OUTPUT."