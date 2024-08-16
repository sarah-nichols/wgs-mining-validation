#!/bin/sh

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load Anaconda3/2024.02-1

# download sequences from wormbase, eupath, ncbi (wgs and nt) and PR2 for building reference databases
set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
# .env DOWNLOAD_SEQS, NCBI_DATASETS_CONDA

#download relevent sequences from ncbi whole genome sequence references
source activate $NCBI_DATASETS_CONDA

set +e
#datasets download genome taxon apicomplexa 
#unzip ncbi_dataset.zip -d "$DOWNLOAD_SEQS"/ncbi/apicomplexa
find "$DOWNLOAD_SEQS"/ncbi/apicomplexa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/apicomplexa \;
cat "$DOWNLOAD_SEQS"/ncbi/apicomplexa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_apicomplexa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/apicomplexa/

datasets download genome taxon cercozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/cercozoa
find "$DOWNLOAD_SEQS"/ncbi/cercozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/cercozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/cercozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_cercozoa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/cercozoa/

datasets download genome taxon euglenozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/euglenozoa
find "$DOWNLOAD_SEQS"/ncbi/euglenozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/euglenozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/euglenozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_euglenozoa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/euglenozoa/

datasets download genome taxon fornicata 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/fornicata
find "$DOWNLOAD_SEQS"/ncbi/fornicata -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/fornicata \;
cat "$DOWNLOAD_SEQS"/ncbi/fornicata/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_fornicata.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/fornicata/

datasets download genome taxon parabasalia 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/parabasalia
find "$DOWNLOAD_SEQS"/ncbi/parabasalia -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/parabasalia \;
cat "$DOWNLOAD_SEQS"/ncbi/parabasalia/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_parabasalia.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/parabasalia/
set -e
conda deactivate