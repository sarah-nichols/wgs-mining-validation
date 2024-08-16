#!/bin/sh

#SBATCH --job-name=random_sub_sample
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 30G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk


module purge
module load SAMtools/1.14-GCC-11.2.0

#set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env
# .env containing HOST_FASTA, HOST_IMG and HOST_HSS

# Desired coverages
coverages=(1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5)

# Loop over each BAM file
for bam_file in "$RAW_DATA"/*.bam; do
  # Get the base name of the BAM file
  bam_base=$(basename "$bam_file" .bam)
  echo "Processing $bam_base"

  # Loop over each desired coverage
  for coverage in "${coverages[@]}"; do
    # Calculate the fraction for subsampling
    echo "Subsampling to $coverage"
    # Create a subsampled BAM file
    samtools view -bs "$coverage" "$bam_file" > "${bam_base}_X${coverage}.bam"
  done
done