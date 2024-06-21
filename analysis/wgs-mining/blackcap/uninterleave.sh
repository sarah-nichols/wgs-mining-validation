#!/bin/sh

# assemble contigs from non-host reads

#SBATCH --job-name=assemble_contigs
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source /data/zool-zost/sann7416/wgs-mining/.env
# .env containing INPUT_FA


for file in ${INT_FASTA}/final.contigs.fa; do

  # Define output file names
  output_r1="${DEINT_FASTA}/$(basename ${file} .fa)_R1.fa"
  output_r2="${DEINT_FASTA}/$(basename ${file} .fa)_R2.fa"

  # Deinterleave the fasta file
  awk 'BEGIN{RS=">"} NR>1 {print ">"$0}' ${file} | paste - - -d '\n' | tee >(cut -f 1-2 -d '>' | tr -d '\n' > ${output_r1}) | cut -f 3- -d '>' | tr -d '\n' > ${output_r2}

done