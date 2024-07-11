#!/bin/sh

#SBATCH --job-name=download_ref
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 10G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load Anaconda3/2024.02-1

#set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

source activate $NCBI_DATASETS_CONDA

datasets download genome accession GCA_009819655.1
unzip ncbi_dataset.zip -d "$HOST_SYLV_PATH"

conda deactivate