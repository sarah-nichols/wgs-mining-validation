#!/bin/sh

#SBATCH --job-name=build_host_ref
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem 80G
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

#set -eu
. "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env" 
# .env containing HOST_FASTA, HOST_IMG and HOST_HSS

echo "Building host reference..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
	-I  "$HOST_SYLV_REG" \
	-O  "$HOST_SYLV_IMG"

echo "Output file: $HOST_SYLV_HSS"

echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$HOST_SYLV_REG" \
   --output "$HOST_SYLV_HSS" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "PathSeqBuildKmers completed"