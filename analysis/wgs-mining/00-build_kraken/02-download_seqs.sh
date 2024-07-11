#!/bin/sh

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

# download sequences from wormbase, eupath and ncbi for building reference genome

set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
# .env DOWNLOAD_SEQS, NCBI_DATASETS_CONDA

#wormbase
wget -P "$DOWNLOAD_SEQS"/wormbase  -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

# Loop through each fasta file, decompress and format header
shopt -s globstar
set +e # Disable exit on error in case some sequences don't work
for file in "$DOWNLOAD_SEQS"/wormbase/**/*.fa.gz; do
  echo "$file"
  gunzip -k "$file"
  decompressed_file="${file%.gz}"
  echo "$decompressed_file"
  # Get the name of the file without the extension
  filename=$(basename "$decompressed_file" .fa)
  echo "$filename"
  bioproject=$(echo ${filename} | grep -oP "PRJ[A-Z]+[0-9]+")
  # Add the file name to the header of each sequence in the file
  sed -i "s/^>/>$bioproject|/" "$decompressed_file"
done
set -e # Re-enable exit on error

#concatanate into a single file
cat $"$DOWNLOAD_SEQS"/wormbase/**/*.fa >  $"$DOWNLOAD_SEQS"/wormbase/wormbase.fasta
rm  $"$DOWNLOAD_SEQS"/wormbase/**/*.fa



#download relevent sequences from eupath
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/AmoebaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/CryptoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/GiardiaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/MicrosporidiaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PiroplasmaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PlasmoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/ToxoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TrichDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TriTrypDB46.tgz

#decompress
for file in "$DOWNLOAD_SEQS"/eupath/*.tgz; do
    tar -xvf "$file" -C "$DOWNLOAD_SEQS"/eupath
done

shopt -s globstar

#concatanate into a single file and delete the rest
cat "$DOWNLOAD_SEQS"/eupath/**/*.fna > "$DOWNLOAD_SEQS"/eupath/eupath.fasta
find "$DOWNLOAD_SEQS"/eupath/ -type f ! -name 'eupath.fasta' -exec rm -f {} +

#format header for kraken2 to recognise
awk '/^>/ {if(seq && ok) {print header ORS seq} header=$0; seq=""; ok = /^>[A-Z0-9]+(\.[0-9]+)? \|/} /^[^>]/ {seq = seq ? seq ORS $0 : $0; ok = ok && !/^N+$/ && !/^x+$/} END {if(seq && ok) print header ORS seq}' $DOWNLOAD_SEQS/eupath/eupath.fasta > $DOWNLOAD_SEQS/eupath/eupath_headers.fasta

shopt -s globstar

#download relevent sequences from ncbi
module load Anaconda3/2024.02-1
source activate $NCBI_DATASETS_CONDA

datasets download genome taxon apicomplexa 
unzip ncbi_dataset.zip -d "$DOWNLOAD_SEQS"/ncbi/apicomplexa
find "$DOWNLOAD_SEQS"/ncbi/apicomplexa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/apicomplexa \;
cat "$DOWNLOAD_SEQS"/ncbi/apicomplexa/*.fna > "$DOWNLOAD_SEQS"/ncbi/apicomplexa.fasta
rm "$DOWNLOAD_SEQS"/ncbi/apicomplexa/*.fna

datasets download genome taxon cercozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/cercozoa
find "$DOWNLOAD_SEQS"/ncbi/cercozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/cercozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/cercozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/cercozoa.fasta
rm "$DOWNLOAD_SEQS"/ncbi/cercozoa/*.fna

datasets download genome taxon euglenozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/euglenozoa
find "$DOWNLOAD_SEQS"/ncbi/euglenozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/euglenozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/euglenozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/euglenozoa.fasta
rm "$DOWNLOAD_SEQS"/ncbi/euglenozoa/*.fna

datasets download genome taxon fornicata 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/fornicata
find "$DOWNLOAD_SEQS"/ncbi/fornicata -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/fornicata \;
cat "$DOWNLOAD_SEQS"/ncbi/fornicata/*.fna > "$DOWNLOAD_SEQS"/ncbi/fornicata.fasta
rm "$DOWNLOAD_SEQS"/ncbi/fornicata/*.fna

datasets download genome taxon parabasalia 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/parabasalia
find "$DOWNLOAD_SEQS"/ncbi/parabasalia -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/parabasalia \;
cat "$DOWNLOAD_SEQS"/ncbi/parabasalia/*.fna > "$DOWNLOAD_SEQS"/ncbi/parabasalia.fasta
rm "$DOWNLOAD_SEQS"/ncbi/parabasalia/*.fna

conda deactivate

wget https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz -O "$DOWNLOAD_SEQS"/pr2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz
gunzip "$DOWNLOAD_SEQS"/pr2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz

