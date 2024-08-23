#!/bin/bash
#SBATCH --job-name=mrun_basta_lca
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --partition medium
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

# Load Anaconda module
module load Anaconda3/2024.02-1

export BASTA_ENV="/data/biol-bird-parasites/sann7416/conda_environments/basta_env"
source activate $BASTA_ENV

#basta taxonomy -d /data/biol-bird-parasites/sann7416/conda_environments/basta_env/.basta/taxonomy

#basta download gb -d /data/biol-bird-parasites/sann7416/conda_environments/basta_env/.basta/taxonomy

basta sequence "$BLAST_OUT/combined_blast_results_fullncbi_100hits_entrezfilter.out" /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/LCA/basta_ncbi_full.tsv gb -d /data/biol-bird-parasites/sann7416/conda_environments/basta_env/.basta/taxonomy -i 97 -p 97 -m 1 -n 0
#basta sequence "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/blast/ncbi_full/combined_blast_results_customncbi.out" /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/LCA/basta_ncbi_full.out gb -d /data/biol-bird-parasites/sann7416/conda_environments/basta_env/.basta/taxonomy -i 95 -p 98 



conda deactivate


