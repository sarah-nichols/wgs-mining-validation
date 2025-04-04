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

#ln -s /data/ncbi/nt/ .

#export BLASTDB=/data/biol-bird-parasites/sann7416/nt


#cd /data/ncbi
# Run BLAST
#blastn -db /data/ncbi/nt -query "$BATCH_FILE" -task blastn -out "$BLAST_OUT/blast_results_full_${SLURM_ARRAY_TASK_ID}.out" -outfmt "6 qseqid sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sskingdoms sscinames"

#blastn -db nt -query "$BATCH_FILE" -task blastn -out "$BLAST_OUT/blast_results_full_${SLURM_ARRAY_TASK_ID}.out" -outfmt "6 qseqid sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" -remote

#blastn -db nt -query "$BATCH_FILE" -task blastn -out "$BLAST_OUT/blast_results_full_${SLURM_ARRAY_TASK_ID}.out" -outfmt "6 qseqid sseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore" -remote

#blastn -query "$BATCH_FILE" -task blastn -entrez_query "all [filter] NOT(environmental samples[All Fields] OR metagenomes[All Fields] OR uncultured[All Fields])" -max_target_seqs 100 -db nt -out "$BLAST_OUT/blast_results_nofilter_100_full_${SLURM_ARRAY_TASK_ID}.out" -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname scomnames sblastname sskingdoms stitle" -remote

blastn -db /data/ncbi/nt/nt -query /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ASVs.fasta -out /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/blast_output_ncbi_full.out
