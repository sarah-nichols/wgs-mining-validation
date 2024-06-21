#!/bin/sh

#SBATCH --job-name=remove_irreg_exp
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env" 
# .env containing HOST_ZOST_IN, HOST_ZOST_OUT

echo "Input file: $HOST_ZOST_IN"
echo "Output file: $HOST_ZOST_OUT"

sed 's/,//' "$HOST_ZOST_IN" > "$HOST_ZOST_OUT"