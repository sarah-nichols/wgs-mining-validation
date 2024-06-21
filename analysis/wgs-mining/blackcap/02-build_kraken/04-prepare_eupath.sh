#!/bin/bash

#SBATCH --job-name=prepare_eupath
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --clusters=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -x
shopt -s globstar
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi
# .env containing EUPATH_PATH

for file in "$EUPATH_PATH"/*.tgz; 
do
  echo "Processing file: $file"
  tar -xzf "$file" -C "$EUPATH_PATH"
done

# Loop through each fasta file
cat $EUPATH_PATH/**/*.fna > eupath.fasta

rm $EUPATH_PATH/**/*.fna

awk '/^>/ {if(seq && ok) {print header ORS seq} header=$0; seq=""; ok = /^>[A-Z0-9]+(\.[0-9]+)? \|/} /^[^>]/ {seq = seq ? seq ORS $0 : $0; ok = ok && !/^N+$/ && !/^x+$/} END {if(seq && ok) print header ORS seq}' $EUPATH_PATH/eupath.fasta > $EUPATH_PATH/eupath_headers.fasta