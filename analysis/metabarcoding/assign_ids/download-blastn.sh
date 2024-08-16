#!/bin/sh

#SBATCH --job-name=blast_download
#SBATCH --clusters=arc
#SBATCH --partition=medium

module load BLAST+/2.14.0-gompi-2022a

export NCBI_nt_DB="/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi/nt"

mkdir $NCBI_nt_DB
cd $NCBI_nt_DB

update_blastdb.pl --decompress nt [*]