#!/bin/sh

#SBATCH --job-name=megablast
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 100000
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk


set -eu
source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

module purge
module load BLAST+/2.14.0-gompi-2022a

#mkdir -p /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ncbi/ncbi_nt
#cd /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ncbi/ncbi_nt
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz

#blastn -db nt -query $ASV_FASTA -taxids 5819 -outfmt 7 -out $BLAST_OUT

blastn -query $ASV_FASTA -db nt -remote -task megablast -entrez_query "haemosporidia[All Fields]" -outfmt 6 -out $BLAST_OUT 

