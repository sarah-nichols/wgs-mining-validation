#!/bin/sh

# download sequences from wormbase, eupath and ncbi for building reference genome

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
# .env DOWNLOAD_SEQS, NCBI_DATASETS_CONDA

#wormbase
wget -P "$DOWNLOAD_SEQS"/wormbase  -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

# Loop through each fasta file
for file in "$DOWNLOAD_SEQS"/wormbase/**/*.fa.gz; 
do
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

cat $"$DOWNLOAD_SEQS"/wormbase/**/*.fa >  $"$DOWNLOAD_SEQS"/wormbase/wormbase.fasta
rm  $"$DOWNLOAD_SEQS"/wormbase/**/*.fa

#eupath
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/AmoebaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/CryptoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/FungiDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/GiardiaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/MicrosporidiaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PiroplasmaDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/PlasmoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/ToxoDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TrichDB46.tgz
wget -P "$DOWNLOAD_SEQS"/eupath ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/TriTrypDB46.tgz

for file in "$DOWNLOAD_SEQS"/eupath/*.tgz; do
    tar -xvf "$file" -C "$DOWNLOAD_SEQS"/eupath
done

cat "$DOWNLOAD_SEQS"/eupath/**/*.fna > "$DOWNLOAD_SEQS"/eupath/eupath.fasta
rm  $"$DOWNLOAD_SEQS"/eupath/**/*.fa

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

conda deactivate $NCBI_DATASETS_CONDA


