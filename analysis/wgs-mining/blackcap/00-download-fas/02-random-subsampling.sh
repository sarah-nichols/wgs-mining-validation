#!/bin/sh

#SBATCH --job-name=random_sub_sample
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 80G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk


module purge
module load SAMtools/1.14-GCC-11.2.0

#set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env
# .env containing HOST_FASTA, HOST_IMG and HOST_HSS

# Directory containing the BAM files
bam_dir="/path/to/bam/files"

# Desired coverages
coverages=(1, 5, 10, 20, 30, 40, 50)

# Loop over each BAM file
for bam_file in "$bam_dir"/*.bam; do
  # Get the base name of the BAM file
  bam_base=$(basename "$bam_file" .bam)

  # Loop over each desired coverage
  for coverage in "${coverages[@]}"; do
    # Calculate the fraction for subsampling
    fraction=$(bc -l <<< "1/$coverage")

    # Create a subsampled BAM file
    samtools view -bs "$fraction" "$bam_file" > "${bam_base}_X${coverage}.bam"
  done
done