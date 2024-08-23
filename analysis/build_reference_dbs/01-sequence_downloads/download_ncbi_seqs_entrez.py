from Bio import Entrez, SeqIO
import time

# Set your email
Entrez.email = "sarahnichols96@outlook.com"

# Define the path to your Taxonomy ID file
taxid_file = "/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/taxids_include.txt"

# Read Taxonomy IDs from the file
with open(taxid_file, "r") as f:
    taxids = [line.strip() for line in f if line.strip()]

# Open a file to write sequences in FASTA format
with open("/data/biol-bird-parasites/sann7416/wgs-mining-validation/data/reference_database/ncbi_custom/ncbi_sequences.fasta", "w") as output_handle:
    for taxid in taxids:
        try:
            # Search for nucleotide sequences by taxonomy ID
            search_handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", retmax=10)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            # Fetch sequences
            if search_results['IdList']:
                fetch_handle = Entrez.efetch(db="nucleotide", id=search_results['IdList'], rettype="fasta", retmode="text")
                sequences = SeqIO.parse(fetch_handle, "fasta")
                SeqIO.write(sequences, output_handle, "fasta")
                fetch_handle.close()
            else:
                print(f"No sequences found for Taxonomy ID {taxid}")
        except Exception as e:
            print(f"Error fetching data for Taxonomy ID {taxid}: {e}")
            time.sleep(5)  # Wait before retrying

print("Sequences have been downloaded and saved to ncbi_sequences.fasta")
