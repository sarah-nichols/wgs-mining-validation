#!/bin/bash
#SBATCH --job-name=combine_blast_results
#SBATCH --output=combine_blast_results.out
#SBATCH --error=combine_blast_results.err
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --time=00:10:00
#SBATCH --mem=1G

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env


# Combine all BLAST results into one file
cat /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/blast/ncbi_full/blast_results_nofilter_*.out > "$BLAST_OUT/combined_blast_results_fullncbi_nofilter_100hits.out"

cat $BLAST_OUT/blast_results_100_full*.out > "$BLAST_OUT/combined_blast_results_fullncbi_100hits_entrezfilter.out"