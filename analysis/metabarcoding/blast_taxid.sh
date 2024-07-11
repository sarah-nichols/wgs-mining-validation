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

blastn -db $NCBI_CUSTOM_DB -query $ASV_FASTA -outfmt 7 -out $BLAST_OUT

#blastn -query $ASV_FASTA -remote -db nt -outfmt 6 -out $BLAST_OUT 

#module load Kraken2/2.1.1-gompi-2020b

#kraken2 --db "$KRAKEN_CUSTOM" \
#  --threads 16 \
#  --classified-out \
#  --use-names \
#  --output "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/metabarcoding.stats.report" \
#  --report "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/metabarcoding.report" \
#   "$ASV_FASTA"


cd "$CONIFER"

# Create an empty file to store the final output
> $CONIFER_OUTPUT_METABARCODING

# Step 1: Extract ASV IDs
grep "^>" $ASV_FASTA | awk '{print substr($0, 2)}' > asv_ids.txt

# Step 2: Run conifer and generate output
./conifer --both_scores -s -i /data/biol-bird-parasites/sann7416/wgs-mining-validation/data/metabarcoding/metabarcoding.stats.report -d $KRAKEN_CUSTOM_TAXID > temp_conifer_output.txt
./conifer --rtl -s -i $KRAKEN_METABARCODING_OUTPUT -d $KRAKEN_CUSTOM_TAXID

# Extract the header from conifer output and modify it to include ASV_ID
head -n 1 temp_conifer_output.txt | sed 's/^/ASV_ID\t/' > header.txt

# Skip the header line from the conifer output for the next steps
tail -n +2 temp_conifer_output.txt > temp_conifer_output_no_header.txt

# Step 3: Combine ASV IDs with conifer output (excluding the original header)
paste asv_ids.txt temp_conifer_output_no_header.txt > combined_output_no_header.txt

# Step 4: Combine the new header with the rest of the output
cat header.txt combined_output_no_header.txt > $CONIFER_OUTPUT_METABARCODING

echo "Number of lines in asv_ids.txt:"
echo "$(wc -l asv_ids.txt)"  # This command counts the number of lines/ASV IDs.

echo "Number of lines in temp_conifer_output.txt:"
echo "$(wc -l temp_conifer_output.txt)"  # Counts the number of lines in the conifer output.

echo "Last 10 lines of asv_ids.txt:"
echo "$(tail asv_ids.txt)"

echo "Last 10 lines of temp_conifer_output_no_header.txt:"
echo "$(tail temp_conifer_output_no_header.txt)"

# Cleanup temporary files
echo "Cleaning up temporary files..."
rm asv_ids.txt temp_conifer_output.txt combined_output_no_header.txt temp_conifer_output_no_header.txt header.txt
echo "Cleanup done."

# Cleanup temporary files
rm asv_ids.txt temp_conifer_output.txt combined_output_no_header.txt temp_conifer_output_no_header.txt header.txt