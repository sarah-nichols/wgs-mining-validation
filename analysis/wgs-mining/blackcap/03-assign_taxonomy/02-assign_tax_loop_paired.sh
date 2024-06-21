#!/bin/sh

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk
#SBATCH --array=1-117%15

set -eu
if [ -f "/data/zool-zost/sann7416/wgs-mining-validation/.env" ]; then
   . "/data/zool-zost/sann7416/wgs-mining-validation/.env"
fi
module purge
module load Kraken2/2.1.1-gompi-2020b
module load SAMtools/1.16.1-GCC-11.3.0
module load SeqKit/2.2.0
module load pear/0.9.11

# Iterate over input BAM files
for bam_file in "$NONHOST_BAMS"/paired_*.bam; do

  # Construct output FASTQ file names
  forward_fastq_file="$NONHOST_PAIRED_F_FASTQS/$(basename "${bam_file%.bam}").1.fastq"
  reverse_fastq_file="$NONHOST_PAIRED_R_FASTQS/$(basename "${bam_file%.bam}").2.fastq"

  # Convert BAM to FASTQ using samtools
  samtools fastq -1 "$forward_fastq_file" -2 "$reverse_fastq_file" -n "$bam_file"

  sample=$(basename "${forward_fastq_file%.1.fastq}.merged" .fastq)

  # Interleave the two FASTQ files into a single FASTQ file using ABYSS
pear -f "$forward_fastq_file" -r "$reverse_fastq_file" -o "${ASSEMBLED_FASTQS}/${sample}".fastq \
 
  kraken2 --db "$KRAKEN_CAPPED_CUSTOM" \
    --threads 16 \
    --classified-out \
    --use-names \
    --output "$KRAKEN_PAIRED_STATS"/"$sample".stats.report \
    --report "$KRAKEN_PAIRED_REPORT"/"$sample".report \
     "${ASSEMBLED_FASTQS}/${sample}"*.assembled.fastq
done