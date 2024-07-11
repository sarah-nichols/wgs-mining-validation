#!/bin/sh

#SBATCH --job-name=build_kraken_db
#SBATCH --partition long
#SBATCH --time=72:00:00
#SBATCH --mem=150G
#SBATCH --cpus-per-task=48  # Request 48 CPUs for this job
#SBATCH --mail-type=ALL
#SBATCH --clusters=all
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env
#env containing SEQS_TO_ADD and KRAKEN_CAPPED paths

module purge
module load Kraken2/2.1.2-gompi-2021b

kraken2-build --download-taxonomy --db "$KRAKEN_CUSTOM2" --threads 48
kraken2-build --download-library protozoa --db "$KRAKEN_CUSTOM2" --threads 48
while IFS= read -r line;
do  
  echo "Adding file: $line to library"
  kraken2-build --add-to-library "$line" --db "$KRAKEN_CUSTOM2" --threads 48 || { echo "kraken2-build failed for $line"; exit 1; }
done < "$SEQS_TO_ADD"

kraken2-build --build --threads 48 --db "$KRAKEN_CUSTOM2" --fast-build
kraken2-build --clean --threads 48 --db "$KRAKEN_CUSTOM2" 




