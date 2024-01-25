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
. "/data/zool-zost/sann7416/wgs-mining-validation/.env" 
# .env containing HOST_FASTA, HOST_IMG and HOST_HSS

echo "Building host reference..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
	-I  "$HOST_FASTA" \
	-O  "$HOST_IMG"

echo "Output file: $HOST_HSS"

echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$HOST_FASTA" \
   --output "$HOST_HSS" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "PathSeqBuildKmers completed"