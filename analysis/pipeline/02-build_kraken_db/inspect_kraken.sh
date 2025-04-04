#!/bin/bash

#SBATCH --job-name=build_kraken_db
#SBATCH --partition short
#SBATCH --mem=20G
#SBATCH --clusters=all
#SBATCH --cpus-per-task=48  # Request 48 CPUs for this job
#SBATCH --mail-type=ALL

module purge
module load Kraken2/2.1.2-gompi-2021b

kraken2-inspect --d /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/data/reference_database/kraken_capped/capped_16gb > /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/data/reference_database/capped16_database_contents.txt