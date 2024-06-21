#!/bin/sh

#SBATCH --job-name=remove_irreg_exp
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source "/data/zool-zost/sann7416/wgs-mining-validation/.env" 
# .env containing HOST_SEQ_IN, HOST_SEQ_OUT

echo "Input file: $HOST_SEQ_IN"
echo "Output file: $HOST_SEQ_OUT"

sed 's/,//' "$HOST_SEQ_IN" > "$HOST_SEQ_OUT"