#!/bin/sh

#SBATCH --job-name=build_kraken_db
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load Kraken2/2.1.2-gompi-2021b

set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"

# .env file contains KRAKEN_CAPPED

wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_16gb_20221209.tar.gz -O "$KRAKEN_CAPPED"