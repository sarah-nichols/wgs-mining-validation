#!/bin/bash

#SBATCH --job-name=build_host_ref
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

#set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env" 
# .env containing HOST_FASTA, HOST_IMG and HOST_HSS

:<<'COMMENT'
echo "Building host reference..."
gatk --java-options "-Xmx60G" BwaMemIndexImageCreator \
	-I  "$HOST_REF_PATH/blackcap_ref_reg.fasta" \
	-O  "$HOST_REF_PATH/blackcap_ref_reg.img"

echo "Output file: $HOST_REF_PATH/blackcap_ref_reg.img"
COMMENT

echo "Running PathSeqBuildKmers..."
gatk --java-options "-Xmx60G" PathSeqBuildKmers  \
   --reference "$HOST_REF_PATH/blackcap_ref_reg.fasta" \
   --output "$HOST_REF_PATH/blackcap_ref_reg.hss" \
   --kmer-mask 16 \
   --kmer-size 31 > pathseq_output.txt 2>&1

echo "PathSeqBuildKmers completed"