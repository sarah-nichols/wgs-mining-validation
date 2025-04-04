#!/bin/sh

#SBATCH --job-name=kraken_assign
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 100000
#SBATCH --clusters=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

module purge
module load Kraken2/2.1.1-gompi-2020b


kraken2 --db "$KRAKEN_CUSTOM" \
  --threads 16 \
  --use-names \
  --output "$KRAKEN_OUTPUT_METABARCODING" \
  --report "$KRAKEN_REPORT_METABARCODING" \
   $ASV_FASTA

cd "$CONIFER"

# Create an empty file to store the final output
> $CONIFER_OUTPUT_METABARCODING
> $CONIFER_OUTPUT_METABARCODING_ASVs


# Print the header line
echo -e "ttaxon_name\ttaxid\treads\tC25\tC50\tC75" > $CONIFER_OUTPUT_METABARCODING

# Run conifer and combine ASV with the output
./conifer --rtl -i "$KRAKEN_OUTPUT_METABARCODING" -d $KRAKEN_CUSTOM_TAXID > $CONIFER_OUTPUT_METABARCODING_ASVs

./conifer -s --rtl -i "$KRAKEN_OUTPUT_METABARCODING" -d $KRAKEN_CUSTOM_TAXID > $CONIFER_OUTPUT_METABARCODING