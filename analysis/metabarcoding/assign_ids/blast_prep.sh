#!/bin/bash

#SBATCH --job-name=blast_batch
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

# Variables
BATCH_SIZE=100
COUNTER=0
BATCH_COUNTER=1

echo "making fastas in $ASV_BATCH"    
# Split FASTA into batches
while IFS= read -r line; do
    if [[ $line == ">"* ]]; then
        COUNTER=$((COUNTER+1))
        if [ $COUNTER -gt $BATCH_SIZE ]; then
            COUNTER=1
            BATCH_COUNTER=$((BATCH_COUNTER+1))
        fi
    fi
    echo "$line" >> "$ASV_BATCH/batch_$BATCH_COUNTER.fasta"
done < "$ASV_FASTA"
