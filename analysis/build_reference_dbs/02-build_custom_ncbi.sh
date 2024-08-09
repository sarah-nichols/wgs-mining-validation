#!/bin/sh

#SBATCH --job-name=build_ncbi_custom_db
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --get-user-env=L
#SBATCH --mem=120G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load BLAST+/2.14.0-gompi-2022a
module load Anaconda3/2024.02-1

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env


# Concatenate fasta files
echo "Concatenating fasta files"

cat "$DOWNLOAD_SEQS/ncbi_custom/ncbi_sequences.fasta" \
    "$DOWNLOAD_SEQS/PR2/pr2_clean.fasta" \
    "$DOWNLOAD_SEQS/eupath/eupath_headers.fasta" > \
    "$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta"
echo "Initial concatenated file sequence count:"
grep -c "^>" "$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta"



# Process fasta file to remove information before and after the accession number in headers and remove duplicates
awk '
BEGIN { seq = ""; header = ""; removed_seqs = 0; removed_dups = 0; }
/^>/ {
    if (seq != "") {
        if (seen[header]++) {
            print header "\n" seq > "'$NCBI_CUSTOM_DB'/removed_duplicates.fasta"
            removed_dups++;
        } else {
            seen[header] = 1;
            print header "\n" seq > "'$NCBI_CUSTOM_DB'/metabarcoding_clean.fasta"
        }
    }
    # Match and extract the accession number
    match($0, /([A-Z]{1,3}\d{5,6}(\.\d+)?|N_?\d{6}(\.\d+)?|NG|NM|NP|NC_\d{6}(\.\d+)?)/, arr);
    if (arr[0] != "") {
        header = ">" arr[0];
        print header > "'$NCBI_CUSTOM_DB'/reformatted_accessions_list.txt"
    } else {
        print $0 > "'$NCBI_CUSTOM_DB'/removed_accessions.fasta"
        removed_seqs++;
    }
    seq = "";
}
!/^>/ {
    seq = seq $0;
}
END {
    if (seq != "") {
        if (seen[header]++) {
            print header "\n" seq > "'$NCBI_CUSTOM_DB'/removed_duplicates.fasta"
            removed_dups++;
        } else {
            print header "\n" seq > "'$NCBI_CUSTOM_DB'/metabarcoding_clean.fasta"
        }
    }
    print "Removed sequences: " removed_seqs > "'$NCBI_CUSTOM_DB'/processing_log.txt"
    print "Removed duplicates: " removed_dups > "'$NCBI_CUSTOM_DB'/processing_log.txt"
}' "$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta"

echo "Cleaned file sequence count:"
grep -c "^>" "$NCBI_CUSTOM_DB/metabarcoding_clean.fasta"
echo "Removed sequences count:"
grep -c "^>" "$NCBI_CUSTOM_DB/removed_accessions.fasta"
echo "Removed duplicates count:"
grep -c "^>" "$NCBI_CUSTOM_DB/removed_duplicates.fasta"

# Download taxid map
#wget -P "$NCBI_CUSTOM_DB"/taxid_map ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#gzip -d "$NCBI_CUSTOM_DB/taxid_map/nucl_gb.accession2taxid.gz"
#sed '1d' "$NCBI_CUSTOM_DB/taxid_map/nucl_gb.accession2taxid" | awk '{print $2" "$3}' > $TAXID_MAP/taxid_map.txt

# Add to custom ncbi db
#cd $NCBI_CUSTOM_DB
#echo "Current directory: $(pwd)"
makeblastdb -in "$NCBI_CUSTOM_DB/metabarcoding_clean.fasta" -title ncbi_parasite_db -dbtype nucl -out "$NCBI_CUSTOM_DB/ncbi_parasite_db" -parse_seqids -taxid_map "$NCBI_CUSTOM_DB/taxid_map.txt"

# Clean up intermediate files
#rm "$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta" "$NCBI_CUSTOM_DB/metabarcoding_clean.fasta"
