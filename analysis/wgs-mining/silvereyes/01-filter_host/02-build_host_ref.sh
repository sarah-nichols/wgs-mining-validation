#!/bin/sh

#SBATCH --job-name=build_host_ref
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 80G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

#set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env" 
# .env containing HOST_ZOST_OUT, HOST_IMG and HOST_ZOST_HSS

echo "Building host reference..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
	-I  "$HOST_ZOST_REG" \
	-O  "$HOST_ZOST_IMG"

echo "Output file: $HOST_ZOST_IMG"

echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$HOST_ZOST_REG" \
   --output "$HOST_ZOST_HSS" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "PathSeqBuildKmers completed"