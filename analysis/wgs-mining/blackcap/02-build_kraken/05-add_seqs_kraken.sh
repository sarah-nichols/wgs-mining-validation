#!/bin/sh

#SBATCH --job-name=add_seqs_kraken_db
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --clusters=all
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi
#env containing SEQS_TO_ADD and KRAKEN_CAPPED paths

module purge
module load Kraken2/2.1.2-gompi-2021b

while IFS= read -r line;
do  
  echo "Adding file: $line to library"
  kraken2-build --add-to-library "$line" --db "$KRAKEN_CAPPED_CUSTOM" || { echo "kraken2-build failed for $line"; exit 1; }
done < "$SEQS_TO_ADD"

kraken2-build --download-taxonomy --db "$KRAKEN_CAPPED_CUSTOM"
kraken2-build --build --db "$KRAKEN_CAPPED_CUSTOM"
kraken2-build --clean --db "$KRAKEN_CAPPED_CUSTOM"




