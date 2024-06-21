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
if [ -f "/data/zool-zost/sann7416/wgs-mining-validation/.env" ]; then
   . "/data/zool-zost/sann7416/wgs-mining-validation/.env"
fi

while IFS= read -r accession
do
    prefetch -O /data/zool-zost/sann7416/wgs-mining-validation/data/wgs-pipeline/coverage_samples "$accession"
done < "$SAMPLE_ACCESSIONS"