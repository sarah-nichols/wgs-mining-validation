#!/bin/bash

#SBATCH --job-name=assemble_contigs_metaspades
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --time=3:00:00
#SBATCH --mem 180G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL

# Script to convert BAM to FASTQ and assemble contigs with metaSPAdes.
# Usage:
#   sbatch --array=1-N%M assembly_metaspades.sh <parameters_file>

set -eu

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <parameters_file>"
    exit 1
fi

PARAMETERS_FILE=$1

# Load parameters from the file
declare -A PARAMS
while IFS="=" read -r key value; do
    if [ -n "$key" ] && [ -n "$value" ]; then
        PARAMS["$key"]="$value"
    fi
done < "$PARAMETERS_FILE"

# Retrieve the required parameters
BAM_PATH=${PARAMS["bam_path"]}
OUTPUT_DIR=${PARAMS["output_dir"]}

# Validate required parameters
if [ -z "$BAM_PATH" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required parameters in the parameters file."
    exit 1
fi

# Load necessary modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load SPAdes/3.15.3-GCC-11.2.0
module load BWA/0.7.17-GCCcore-11.2.0
module load seqtk/1.3-GCC-11.2.0

# Create output directories
PAIRED_F_FASTQS="$OUTPUT_DIR/paired_fastqs_forward"
PAIRED_R_FASTQS="$OUTPUT_DIR/paired_fastqs_reverse"
SPADES_OUTPUT="$OUTPUT_DIR/spades_output"
UNASSEMBLED_READS="$OUTPUT_DIR/unassembled_reads"

mkdir -p "$PAIRED_F_FASTQS"
mkdir -p "$PAIRED_R_FASTQS"
mkdir -p "$SPADES_OUTPUT"
mkdir -p "$UNASSEMBLED_READS"

# List all BAM files
BAM_FILES=("$BAM_PATH"/paired_*.bam)
NUM_FILES=${#BAM_FILES[@]}

if [ "$NUM_FILES" -eq 0 ]; then
    echo "Error: No BAM files found."
    exit 1
fi

# Use SLURM_ARRAY_TASK_ID to pick the specific file
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
bam_file=${BAM_FILES[$INDEX]}

if [ -z "$bam_file" ]; then
    echo "Error: No BAM file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
    exit 1
fi

# Generate FASTQ filenames
forward_fastq_file="$PAIRED_F_FASTQS/$(basename "${bam_file%.bam}").1.fastq"
reverse_fastq_file="$PAIRED_R_FASTQS/$(basename "${bam_file%.bam}").2.fastq"

# Convert BAM to FASTQ using samtools
samtools fastq -1 "$forward_fastq_file" -2 "$reverse_fastq_file" -n "$bam_file"

# Assemble contigs using metaSPAdes
spades_output_path="$SPADES_OUTPUT/$(basename "${bam_file%.bam}")"
spades.py --meta -1 "$forward_fastq_file" -2 "$reverse_fastq_file" -o "$spades_output_path"

# Step 1: Index the contigs for BWA-MEM
contigs_file="$spades_output_path/contigs.fasta"
bwa index "$contigs_file"

# Step 2: Align reads to contigs using BWA-MEM
bwa mem -t 16 "$contigs_file" "$forward_fastq_file" "$reverse_fastq_file" > "$OUTPUT_DIR/$(basename "${bam_file%.bam}").sam"

# Step 3: Convert SAM to BAM and extract unmapped reads
samtools view -bS -F 4 "$OUTPUT_DIR/$(basename "${bam_file%.bam}").sam" > "$OUTPUT_DIR/$(basename "${bam_file%.bam}").mapped.bam"
samtools view -bS -f 4 "$OUTPUT_DIR/$(basename "${bam_file%.bam}").sam" > "$OUTPUT_DIR/$(basename "${bam_file%.bam}").unmapped.bam"

# Step 4: Convert unmapped BAM to FASTQ
samtools fastq -1 "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.1.fastq" \
               -2 "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.2.fastq" \
               "$OUTPUT_DIR/$(basename "${bam_file%.bam}").unmapped.bam"

# Step 5: Convert unassembled FASTQ to FASTA
seqtk seq -A "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.1.fastq" > "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.1.fasta"
seqtk seq -A "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.2.fastq" > "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.2.fasta"

# Step 6: Combine contigs and unassembled reads into one FASTA file
combined_fasta="$OUTPUT_DIR/combined_sequences/$(basename "${bam_file%.bam}").combined_sequences.fasta"
mkdir -p "$(dirname "$combined_fasta")"
cat "$contigs_file" "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.1.fasta" "$UNASSEMBLED_READS/$(basename "${bam_file%.bam}").unassembled.2.fasta" > "$combined_fasta"
