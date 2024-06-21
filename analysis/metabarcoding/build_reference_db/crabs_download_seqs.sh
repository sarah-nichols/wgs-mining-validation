#!/bin/bash

#SBATCH --job-name=download_ref_seqs
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu

module load Anaconda3/2024.02-1

# Activate Conda environment
export CONPREFIX='/data/biol-bird-parasites/sann7416/conda_environments/CRABS_python3.8'
conda activate $CONPREFIX

crabs db_download --source ncbi --database nucleotide --query '((haemosporidian[All Fields] NOT cytB[All Fields]) NOT cytochrome[All Fields]) NOT cox[All Fields] NOT ("Chordata"[Organism] OR Chordata[All Fields]) NOT "Nycteribia schmidlii"[porgn]' --output /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ncbi/haemosporidians_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
crabs db_download --source ncbi --database gene --query '"Haemosporida"[Organism] AND 18S[All Fields]' --output /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ncbi/haemosporidians_gene.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
crabs db_download --source ncbi --database genome --query '"Haemosporida"[Organism]' --output /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ncbi/haemosporidians_genome.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000

