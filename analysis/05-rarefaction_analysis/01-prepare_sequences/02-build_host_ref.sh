#!/bin/bash

#SBATCH --job-name=build_host_ref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition short
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL

# Script to build host reference using GATK's BwaMemIndexImageCreator and PathSeqBuildKmers, and optionally run bwa index.
# Usage:
#   sbatch this_script.sh <parameters_file> [--mail-user=email@example.com]
#
# Arguments:
#   <parameters_file> - Path to the text file containing paths to the host reference directory
#   [email]           - (Optional) Email address to receive job notifications

module purge
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8
module load BWA/0.7.17-GCCcore-11.2.0

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
HOST_REF_DIR=${PARAMS["host_ref_dir"]}

# Validate that the required parameter is provided
if [ -z "$HOST_REF_DIR" ]; then
    echo "Error: Missing required parameter for host reference directory in the parameters file."
    exit 1
fi

# Index the reference genome using bwa index
echo "Indexing the reference genome using bwa index..."
bwa index "$HOST_REF_DIR/blackcap_ref_reg.fasta"

#echo "BWA index completed"

# Build host reference index using BwaMemIndexImageCreator
echo "Building host reference index (IMG)..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
    -I "$HOST_REF_DIR/blackcap_ref_reg.fasta" \
    -O "$HOST_REF_DIR/blackcap_ref_reg.img"

echo "Output file: $HOST_REF_DIR/blackcap_ref_reg.img"

# Build host reference using PathSeqBuildKmers
echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$HOST_REF_DIR/blackcap_ref_reg.fasta" \
   --output "$HOST_REF_DIR/blackcap_ref_reg.hss" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "Output file: $HOST_REF_DIR/blackcap_ref_reg.hss"
echo "PathSeqBuildKmers completed"
