#!/bin/bash

#SBATCH --job-name=download_samples
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 64000
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load SRA-Toolkit/3.0.3-gompi-2022a

set -x
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

while IFS= read -r accession
do
    prefetch -O $RAW_DATA "$accession"
    sam-dump "$accession" | samtools view -bS - > "${RAW_DATA}/${accession}.bam"
done < "$SAMPLE_ACCESSIONS_SYLV"

