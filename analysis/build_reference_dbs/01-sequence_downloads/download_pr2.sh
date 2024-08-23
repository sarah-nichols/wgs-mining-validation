#!/bin/sh

#SBATCH --job-name=download_pr2
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

# download sequences from wormbase, eupath, ncbi (wgs and nt) and PR2 for building reference databases
set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"

#download PR2 database
wget -O "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
gzip -d "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" > "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta"

#filter sequences from PR2 that are not of interest
awk 'BEGIN{IGNORECASE=1} /^>/{printit=1; for (i=1;i<=NF;i++) {if ($i ~ /16S_rRNA|Fungi|uncultured|Arthropoda|Annelida|Craniata|Streptophyta|Chlorophyta/) {printit=0}} } {if (printit) print}' "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta" > "$DOWNLOAD_SEQS/PR2/pr2_clean.fasta"

