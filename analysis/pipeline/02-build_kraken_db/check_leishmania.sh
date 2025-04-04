#!/bin/bash

#SBATCH --job-name=leishmania_extract
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 200G
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL

module load Kraken2/2.1.1-gompi-2020b

kraken2-inspect --db /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/data/reference_database/kraken_custom | grep -w "5666" > /data/biol-bird-parasites/sann7416/island-biogeography-wgs-mining/data/reference_database/taxon_5666_info.txt