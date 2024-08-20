#!/bin/sh

#SBATCH --job-name=download_wormbase
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=10G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

# download sequences from wormbase, eupath, ncbi (wgs and nt) and PR2 for building reference databases
set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
# .env DOWNLOAD_SEQS, NCBI_DATASETS_CONDA

#wormbase
wget -P "$DOWNLOAD_SEQS"/wormbase  -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

# Loop through each fasta file, decompress and format header
shopt -s globstar
set +e # Disable exit on error in case some sequences don't work
for file in "$DOWNLOAD_SEQS"/wormbase/**/*.fa.gz; do
  echo "$file"
  gunzip -k "$file"
  decompressed_file="${file%.gz}"
  echo "$decompressed_file"
  # Get the name of the file without the extension
  filename=$(basename "$decompressed_file" .fa)
  echo "$filename"
  bioproject=$(echo ${filename} | grep -oP "PRJ[A-Z]+[0-9]+")
  # Add the file name to the header of each sequence in the file
  sed -i "s/^>/>$bioproject|/" "$decompressed_file"
done
set -e # Re-enable exit on error

#concatanate into a single file
cat $"$DOWNLOAD_SEQS"/wormbase/**/*.fa >  $"$DOWNLOAD_SEQS"/wormbase/wormbase.fasta

rm -r "$DOWNLOAD_SEQS/wormbase/ftp.ebi.ac.uk"