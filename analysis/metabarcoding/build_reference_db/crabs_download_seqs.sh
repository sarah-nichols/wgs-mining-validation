#!/bin/bash

#SBATCH --job-name=download_ref_seqs
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --mem 64000
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu

module load Anaconda3/2024.02-1

# Activate Conda environment
export CONPREFIX='/data/zool-zost/sann7416/crabenv'
source activate $CONPREFIX

cd /data/zool-zost/sann7416/reference_database_creator-main/

./crabs db_download --source ncbi --database nucleotide --query '((haemosporidian[All Fields] NOT cytB[All Fields]) NOT cytochrome[All Fields]) NOT cox[All Fields] NOT ("Chordata"[Organism] OR Chordata[All Fields]) NOT "Nycteribia schmidlii"[porgn]' --output haemosporidians_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000

./crabs db_download --source taxonomy

./crabs insilico_pcr --input haemosporidians_nucleotide.fasta --output haemosporidians_nucleotide_primers.fasta --fwd GCAAGTCTGGTGCCAG --rev CTTTAARTTTCASYCTTGCG --error 4.5

./crabs assign_tax --input haemosporidians_nucleotide_primers.fasta --output haemosporidians_nucleotide_primers.tsv --acc2tax nucl_gb.accession2taxid --taxid nodes.dmp --name names.dmp --missing missing_taxa.tsv

./crabs visualization --method diversity --input haemosporidians_nucleotide_primers.tsv --level species --output haemosporidians_nucleotide_primers_species_diversity.png

./crabs tax_format --input haemosporidians_nucleotide_primers.tsv --output haemosporidian_references.fsa --format sintax
./crabs tax_format --input haemosporidians_nucleotide_primers.tsv --output haemosporidian_references.fasta.gz --format idt

cp haemosporidian_references.fasta.gz /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/reference_database/haemosporidian_reference.fasta.gz
vsearch --sintax /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/ASVs.fasta --db haemosporidian_references.fsa --tabbedout vsearch_output.txt 