#!/bin/sh

#SBATCH --job-name=build_ncbi_custom_db
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --get-user-env=L
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load BLAST+/2.14.0-gompi-2022a
module load Anaconda3/2024.02-1

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env
#CONTAINS NCBI_NT_PARASITES, PR2, EUPATH, PARASITE_FASTA, NCBI_CUSTOM_DB

# Activate Conda environment
export CONPREFIX='/data/zool-zost/sann7416/crabenv'
source activate $CONPREFIX

cd /data/zool-zost/sann7416/reference_database_creator-main/

./crabs db_download --source ncbi --database nucleotide --query '(((((((("Apicomplexa"[Organism] OR Apicomplexa[All Fields]) OR ("Parabasalia"[Organism] OR Parabasalia[All Fields])) OR ("Fornicata"[Organism] OR Fornicata[All Fields])) OR ("Euglenozoa"[Organism] OR Euglenozoa[All Fields])) OR ("Cercozoa"[Organism] OR Cercozoa[All Fields])) AND Ribosomal[All Fields]) NOT uncultured[All Fields]) NOT environmental[All Fields]) NOT ("unclassified sequences"[Organism] OR unclassified[All Fields]) AND biomol_rrna[PROP]' --output protist_parasite_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
cp protist_parasite_nucleotide.fasta $NCBI_NT_PARASITES

conda deactivate

#concatanate into one fasta
echo "Concatenating fasta files"
cat $NCBI_NT_PARASITES $PR2 $EUPATH > $PARASITE_FASTA_REPEATS

#remove information after the accession
awk '/^>/ {
    # Check if there is a period in the header
    pos = index($0, ".");
    if (pos > 0 && length($0) > pos) {
        # If a period exists, keep up to and including the first character after the period
        $0 = substr($0, 1, pos+1);
    }
    # No else block needed; if no period, the header remains unchanged
    if (seq != "") print header"\n"seq;
    header=$0; seq="";
}
!/^>/ {
    seq = seq $0;
}
END {
    if (seq != "") print header"\n"seq;
}' $PARASITE_FASTA_REPEATS > $PARASITE_FASTA_FORMAT

#remove duplicates
awk '/^>/{key=tolower($0); if(seen[key]++) skip=1; else {print; skip=0}} !/^>/{if(!skip) print}' $PARASITE_FASTA_FORMAT > $PARASITE_FASTA


# Add to custom ncbi db
makeblastdb -in $PARASITE_FASTA -title ncbi_parasite_db -dbtype nucl -out $NCBI_CUSTOM_DB -parse_seqids