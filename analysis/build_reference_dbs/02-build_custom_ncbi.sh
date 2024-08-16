#!/bin/sh

#SBATCH --job-name=build_ncbi_custom_db
#SBATCH --nodes=1
#SBATCH --partition short
#SBATCH --get-user-env=L
#SBATCH --mem=30G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sarah.nichols@biology.ox.ac.uk

module purge
module load BLAST+/2.14.0-gompi-2022a
module load Anaconda3/2024.02-1
module load SeqKit/2.2.0

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

# Define file paths
CONCATENATED_FASTA="$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta"
REFORMATTED_FASTA="$NCBI_CUSTOM_DB/metabarcoding_accessions_reformatted.fasta"
REFORMATTED_ACCESSIONS_LIST="$NCBI_CUSTOM_DB/reformatted_accessions_list.txt"
REMOVED_DUPLICATES_FASTA="$NCBI_CUSTOM_DB/metabarcoding_duplicates.fasta"
CLEAN_FASTA="$NCBI_CUSTOM_DB/metabarcoding_clean.fasta"

# Concatenate fasta files
echo "Concatenating fasta files"

# Concatenate the three FASTA files into a new file
cat "$DOWNLOAD_SEQS/ncbi_custom/ncbi_sequences.fasta" \
    "$DOWNLOAD_SEQS/PR2/pr2_clean.fasta" \
    "$DOWNLOAD_SEQS/eupath/eupath_headers.fasta" > "$CONCATENATED_FASTA"

# Print the number of headers in the concatenated file
echo "Initial concatenated file sequence count:"
grep -c "^>" "$CONCATENATED_FASTA"

# Reformat headers to retain only the NCBI accession numbers and discard sequences without accession numbers
awk '
BEGIN {
    print "Reformatting headers to retain only accession numbers and discarding sequences without accession numbers.";
    header_count = 0;
    discarded_count = 0;
    seq = "";
}
# Process header lines
/^>/ {
    # If a previous sequence was collected, check the header for accession and process it
    if (seq != "" && header != "") {
        print header "\n" seq >> "'$REFORMATTED_FASTA'";
        seq = "";  # Reset sequence buffer
    }

    # Extract the accession number from the header
    # The regex below matches typical NCBI accession numbers
    # Example formats include: ABC12345, ABC12345.1, XP_123456, XP_123456.1, etc.
    if (match($0, /([A-Z]{1,3}_?[0-9]+\.[0-9]+|[A-Z]{1,3}[0-9]{5,6}(\.[0-9]+)?)/, arr)) {
        header = ">" arr[1];
        print arr[1] >> "'$REFORMATTED_ACCESSIONS_LIST'"
        header_count++;
    } else {
        # If no accession found, discard this sequence
        header = "";  # Clear the header as a flag for discard
        discarded_count++;
    }
    next;
}

# Process sequence lines
{
    if (header != "") {
        seq = seq $0;
    }
}

END {
    # Check the last sequence and write it if the header is valid
    if (seq != "" && header != "") {
        print header "\n" seq >> "'$REFORMATTED_FASTA'";
    }
    print "Total reformatted headers: " header_count > "/dev/stderr";
    print "Total discarded sequences: " discarded_count > "/dev/stderr";
}
' "$CONCATENATED_FASTA"

echo "Reformatted headers have been saved to: $REFORMATTED_FASTA"
echo "List of reformatted accessions has been saved to: $REFORMATTED_ACCESSIONS_LIST"

# Remove duplicate sequences
seqkit rmdup -n "$REFORMATTED_FASTA" -o "$CLEAN_FASTA" -D "$REMOVED_DUPLICATES_FASTA"

# Download taxid map
#wget -P "$NCBI_CUSTOM_DB"/taxid_map ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
#gzip -d "$NCBI_CUSTOM_DB/taxid_map/nucl_gb.accession2taxid.gz"
#sed '1d' "$NCBI_CUSTOM_DB/taxid_map/nucl_gb.accession2taxid" | awk '{print $2" "$3}' > $TAXID_MAP/taxid_map.txt

# Add to custom ncbi db
cd $NCBI_CUSTOM_DB
#echo "Current directory: $(pwd)"
makeblastdb -in "$CLEAN_FASTA" -title ncbi_parasite_db -dbtype nucl -out "$NCBI_CUSTOM_DB/ncbi_parasite_db" -parse_seqids -taxid_map "$NCBI_CUSTOM_DB/taxid_map.txt"

# Clean up intermediate files
#rm "$NCBI_CUSTOM_DB/metabarcoding_accessions.fasta" "$NCBI_CUSTOM_DB/metabarcoding_clean.fasta"
