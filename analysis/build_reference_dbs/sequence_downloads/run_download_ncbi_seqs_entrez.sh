#!/bin/bash
#SBATCH --job-name=my_python_job
#SBATCH --output=output.txt
#SBATCH --error=error.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition long
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

# Load Anaconda module
module load Anaconda3/2024.02-1

TAXKIT_ENV="/data/biol-bird-parasites/sann7416/conda_environments/taxonkit_env"
source activate $TAXKIT_ENV

taxonkit list \
    --ids 5794,207245,5719,33682,136419 \
    --indent "" --data-dir /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/taxid_map/taxdump \
    --show-name |
     grep -v -E 'uncultured|environmental' | 
     awk '{print $1}' > /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/taxids_include.txt

conda deactivate 

export CONPREFIX=/data/biol-bird-parasites/sann7416/conda_environments/entrez_env

# Activate your conda environment
source activate $CONPREFIX

# Run the Python script
python /data/biol-bird-parasites/sann7416/wgs-mining-validation/analysis/build_reference_dbs/sequence_downloads/download_ncbi_seqs_entrez.py

# Deactivate the environment
conda deactivate
