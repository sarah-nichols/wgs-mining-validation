#!/bin/bash

#SBATCH --job-name=prepare_wormbase
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -x
shopt -s globstar
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi
# .env containing WORMBASE_PATH
echo "$WORMBASE_PATH"

# Loop through each fasta file
for file in "$WORMBASE_PATH"/**/*.fa.gz; 
do
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

cat $WORMBASE_PATH/**/*.fa > $WORMBASE_PATH/wormbase.fasta

rm "$WORMBASE_PATH"/**/*.fa
