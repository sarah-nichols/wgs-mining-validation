#!/bin/sh

#SBATCH --job-name=download_kraken_capped
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load Kraken2/2.1.2-gompi-2021b

set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"

# .env file contains KRAKEN_CAPPED
echo "Downloading to: $KRAKEN_CAPPED"
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz -O "$KRAKEN_CAPPED"

tar -xvf "$KRAKEN_CAPPED" -C "$KRAKEN_CUSTOM"