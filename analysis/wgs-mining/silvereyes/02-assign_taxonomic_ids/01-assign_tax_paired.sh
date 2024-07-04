#!/bin/sh

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk
#SBATCH --array=1-72%10

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

module purge
module load Kraken2/2.1.1-gompi-2020b
module load SAMtools/1.16.1-GCC-11.3.0
module load SeqKit/2.2.0
module load pear/0.9.11

# List all BAM files in the directory
BAM_FILES=("$BAM_ZOST_PATH"/paired/*.bam)

# Use SLURM_ARRAY_TASK_ID to pick a specific file
# Note: Bash arrays are 0-indexed, but SLURM_ARRAY_TASK_ID starts at 1
INDEX=$((SLURM_ARRAY_TASK_ID - 1))
bam_file=${BAM_FILES[$INDEX]}

# Ensure bam_file is not empty or unset
if [ -z "$bam_file" ]; then
  echo "Error: No BAM file found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
  exit 1
fi

# Construct output FASTQ file names
forward_fastq_file="$PAIRED_F_FASTQS_ZOST/$(basename "${bam_file%.bam}").1.fastq"
reverse_fastq_file="$PAIRED_R_FASTQS_ZOST/$(basename "${bam_file%.bam}").2.fastq"

# Convert BAM to FASTQ using samtools
samtools fastq -1 "$forward_fastq_file" -2 "$reverse_fastq_file" -n "$bam_file"

sample=$(basename "${forward_fastq_file%.1.fastq}.merged" .fastq)

# Interleave the two FASTQ files into a single FASTQ file using ABYSS
pear -f "$forward_fastq_file" -r "$reverse_fastq_file" -o "${ASSEMBLED_FASTQS_ZOST}/${sample}".fastq

 
kraken2 --db "$KRAKEN_CUSTOM" \
  --threads 16 \
  --classified-out \
  --use-names \
  --output "$KRAKEN_OUTPUT_ZOST_PAIRED"/"$sample".stats.report \
  --report "$KRAKEN_REPORT_ZOST_PAIRED"/"$sample".report \
   "${ASSEMBLED_FASTQS_ZOST}/${sample}"*.assembled.fastq