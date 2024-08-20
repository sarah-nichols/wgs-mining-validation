#!/bin/sh

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=120G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

#wget -P "$DOWNLOAD_SEQS"/megan https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-nucl-Feb2022.db.zip
unzip "$DOWNLOAD_SEQS"/megan/megan-nucl-Feb2022.db.zip -d "$DOWNLOAD_SEQS"/megan
