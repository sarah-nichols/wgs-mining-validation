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
unzip ncbi_dataset.zip -d "$HOST_PATH/blackcap_GCA_009819655.1"

conda deactivate

echo "Input file: "$HOST_REF_PATH/blackcap_GCA_009819655.1/ncbi_dataset/data/GCA_009819655.1/GCA_009819655.1_bSylAtr1.pri_genomic.fna""
echo "Output file: "$HOST_REF_PATH/blackcap_ref_reg.fasta""

sed 's/,//' "$HOST_REF_PATH/blackcap_GCA_009819655.1/ncbi_dataset/data/GCA_009819655.1/GCA_009819655.1_bSylAtr1.pri_genomic.fna" > "$HOST_PATH/blackcap_ref_reg.fasta"

rm -r "$HOST_REF_PATH/blackcap_GCA_009819655.1"