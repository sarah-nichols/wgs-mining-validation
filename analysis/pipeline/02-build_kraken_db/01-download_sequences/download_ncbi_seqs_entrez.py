import argparse
from Bio import Entrez, SeqIO
import time

# Define a function to parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Download sequences using NCBI Entrez API based on Taxonomy IDs")
    parser.add_argument('--taxid-file', required=True, help="Path to the file containing Taxonomy IDs")
    parser.add_argument('--output', required=True, help="Path to the output FASTA file")
    args = parser.parse_args()
    return args

# Function to download sequences from NCBI
def download_sequences(taxids, output_path):
    # Set the Entrez email (important for NCBI API usage)
    Entrez.email = "sarahnichols96@outlook.com"
    
    # Open output file for writing sequences
    with open(output_path, "w") as output_handle:
        # Loop through the list of Taxonomy IDs and fetch sequences
        for taxid in taxids:
            try:
                # Search for nucleotide sequences by taxonomy ID
                search_handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", retmax=10)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                # Fetch the sequences using the IDs returned from the search
                if search_results['IdList']:
                    fetch_handle = Entrez.efetch(db="nucleotide", id=search_results['IdList'], rettype="fasta", retmode="text")
                    sequences = SeqIO.parse(fetch_handle, "fasta")
                    SeqIO.write(sequences, output_handle, "fasta")
                    fetch_handle.close()
                else:
                    print(f"No sequences found for Taxonomy ID {taxid}")
            except Exception as e:
                print(f"Error fetching data for Taxonomy ID {taxid}: {e}")
                time.sleep(5)  # Wait 5 seconds before retrying in case of API issues

    print(f"Sequences have been downloaded and saved to {output_path}")

# Main execution flow
if __name__ == "__main__":
    # Parse the command-line arguments
    args = parse_args()

    # Read the Taxonomy IDs from the provided file
    with open(args.taxid_file, "r") as f:
        taxids = [line.strip() for line in f if line.strip()]

    # Download sequences and write them to the output FASTA file
    download_sequences(taxids, args.output)
