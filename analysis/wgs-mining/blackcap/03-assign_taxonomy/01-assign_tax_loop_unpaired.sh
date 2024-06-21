#!/bin/sh

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi
module purge
module load Kraken2/2.1.1-gompi-2020b
module load SAMtools/1.16.1-GCC-11.3.0

# Iterate over input BAM files
for bam_file in "$NONHOST_BAMS"/unpaired_*.bam; do

  # Construct output FASTA file name
  fasta_file="$NONHOST_UNPAIRED_FASTAS/$(basename "${bam_file%.bam}").fasta"
  # Convert BAM to FASTA using samtools
  samtools fasta "$bam_file" > "$fasta_file"
  
  sample=$(basename "$fasta_file" .fasta)
  kraken2 --db "$KRAKEN_CAPPED_CUSTOM" \
    --threads 16 \
    --classified-out \
    --use-names \
    --output "$KRAKEN_UNPAIRED_STATS"/"$sample".stats.report \
    --report "$KRAKEN_UNPAIRED_REPORT"/"$sample".report \
    "$fasta_file"

done


