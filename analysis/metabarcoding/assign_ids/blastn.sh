#!/bin/sh

#SBATCH --job-name=blast_search
#SBATCH --output=/data/biol-bird-parasites/sann7416/blast_%A_%a.out
#SBATCH --error=/data/biol-bird-parasites/sann7416/blast_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --clusters=all
#SBATCH --partition=short
#SBATCH --array=1-11  # N is the total number of batches

# Load BLAST module
module load BLAST+/2.14.0-gompi-2022a
module load Anaconda3/2024.02-1

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

# Variables
BATCH_FILE="$ASV_BATCH/batch_${SLURM_ARRAY_TASK_ID}.fasta"
echo "BLASTing $BATCH_FILE"

cd $NCBI_CUSTOM_DB
# Run BLAST
blastn -db ncbi_parasite_db -query "$BATCH_FILE" -task blastn -out "$BLAST_OUT/blast_results_${SLURM_ARRAY_TASK_ID}.out" -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomnames sblastname sskingdoms stitle"

#blastn -db ncbi_parasite_db -query /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ASVs.fasta -task blastn -out "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/blast/all_blast.out.tab" -outfmt "6 qseqid sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sskingdoms sscinames"
