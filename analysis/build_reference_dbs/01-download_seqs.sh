#!/bin/sh

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

module purge
module load Anaconda3/2024.02-1

# download sequences from wormbase, eupath, ncbi (wgs and nt) and PR2 for building reference databases
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

rm -r "$DOWNLOAD_SEQS/wormbase/ftp.ebi.ac.uk"

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

#format header for kraken2/ ncbi to recognise
awk '/^>/ {if(seq && ok) {print header ORS seq} header=$0; seq=""; ok = /^>[A-Z0-9]+(\.[0-9]+)? \|/} /^[^>]/ {seq = seq ? seq ORS $0 : $0; ok = ok && !/^N+$/ && !/^x+$/} END {if(seq && ok) print header ORS seq}' $DOWNLOAD_SEQS/eupath/eupath.fasta > $DOWNLOAD_SEQS/eupath/eupath_headers.fasta

rm -rf "$DOWNLOAD_SEQS"/eupath/AmoebaDB "$DOWNLOAD_SEQS"/eupath/CryptoDB "$DOWNLOAD_SEQS"/eupath/GiardiaDB "$DOWNLOAD_SEQS"/eupath/MicrosporidiaDB "$DOWNLOAD_SEQS"/eupath/PiroplasmaDB "$DOWNLOAD_SEQS"/eupath/PlasmoDB "$DOWNLOAD_SEQS"/eupath/ToxoDB "$DOWNLOAD_SEQS"/eupath/TrichDB "$DOWNLOAD_SEQS"/eupath/TriTrypDB

shopt -s globstar

# Activate Conda environment

:<<'COMMENT'
source activate $CRAB_ENV

cd /data/zool-zost/sann7416/reference_database_creator-main/

#Download sequences that match entrez query
./crabs db_download --source ncbi --database nucleotide --query 'Haemoproteus "Haemoproteus"[Organism] OR haemoproteus[All Fields] AND ribosom*' --output haemoproteus_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query '"Leucocytozoon"[Organism] OR Leucocytozoon[All Fields] AND ribosomal' --output leucocytozoon_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Plasmodium[All Fields] AND ribosomal' --output plasmodium_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Atoxoplasma[All Fields] AND ribosomal' --output atoxoplasma_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Isospora[All Fields] AND ribosomal' --output isospora_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Giardia[All Fields] AND ribosomal' --output giardia_parasite_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Trichomonas[All Fields] AND ribosomal' --output trichomonas_parasite_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Hexamita[All Fields] OR Spironucleus[All Fields] OR Octomitus [All Fields] AND ribosomal' --output hexamita_parasite_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Eimeria[All Fields] AND ribosomal NOT uncultured bacterium' --output eimeria_parasite_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Sarcocystis[All Fields]  AND ribosomal' --output sarcocystis_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'cryptosporidium[All Fields] AND ribosomal NOT environmental' --output cryptosporidium_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Sternostoma[All Fields] AND ribosomal' --output sternostoma_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'trypanosoma[All Fields] AND ribosomal' --output trypanosoma_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Babesia[All Fields] AND ribosomal' --output babesia_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Theileria[All Fields] AND ribosomal' --output theileria_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Hepatozoon[All Fields] AND ribosomal' --output hepatozoon_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Cytauxzoon[All Fields] AND ribosomal' --output cytauxzoon_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Besnoitia[All Fields] AND ribosomal' --output besnoitia_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Neospora[All Fields]' --output neospora_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000
./crabs db_download --source ncbi --database nucleotide --query 'Toxoplasma AND ribosomal' --output toxoplasma_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000

./crabs db_download --source ncbi --database nucleotide --query  --output haemoproteus_nucleotide.fasta --keep_original yes --email sarahnichols96@outlook.com --batchsize 5000


cat haemoproteus_nucleotide.fasta leucocytozoon_nucleotide.fasta plasmodium_nucleotide.fasta atoxoplasma_nucleotide.fasta isospora_nucleotide.fasta giardia_parasite_nucleotide.fasta trichomonas_parasite_nucleotide.fasta hexamita_parasite_nucleotide.fasta eimeria_parasite_nucleotide.fasta sarcocystis_nucleotide.fasta cryptosporidium_nucleotide.fasta sternostoma_nucleotide.fasta trypanosoma_nucleotide.fasta babesia_nucleotide.fasta theileria_nucleotide.fasta hepatozoon_nucleotide.fasta cytauxzoon_nucleotide.fasta besnoitia_nucleotide.fasta neospora_nucleotide.fasta toxoplasma_nucleotide.fasta > protist_parasite_nucleotide.fasta
cp protist_parasite_nucleotide.fasta $DOWNLOAD_SEQS/ncbi/protist_parasite_nucleotide.fasta

conda deactivate
COMMENT



echo "export PATH=/data/biol-bird-parasites/sann7416/wgs-mining-validation/software/edirect:\$PATH" >> ~/.bashrc
export PATH=${EDIRECT_PATH}/edirect:${PATH}

export PATH=$PATH:/data/biol-bird-parasites/sann7416/wgs-mining-validation/software/edirect



while read taxid; do
    esearch -db nucleotide -query "txid${taxid}[Organism]" -email sarahnichols96@outlook.com | efetch -format fasta >> ncbi_sequences.fasta
done < "$TAXID_LIST"


#download relevent sequences from ncbi whole genome sequence references
source activate $NCBI_DATASETS_CONDA

set +e
#datasets download genome taxon apicomplexa 
#unzip ncbi_dataset.zip -d "$DOWNLOAD_SEQS"/ncbi/apicomplexa
find "$DOWNLOAD_SEQS"/ncbi/apicomplexa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/apicomplexa \;
cat "$DOWNLOAD_SEQS"/ncbi/apicomplexa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_apicomplexa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/apicomplexa/

datasets download genome taxon cercozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/cercozoa
find "$DOWNLOAD_SEQS"/ncbi/cercozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/cercozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/cercozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_cercozoa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/cercozoa/

datasets download genome taxon euglenozoa
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/euglenozoa
find "$DOWNLOAD_SEQS"/ncbi/euglenozoa -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/euglenozoa \;
cat "$DOWNLOAD_SEQS"/ncbi/euglenozoa/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_euglenozoa.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/euglenozoa/

datasets download genome taxon fornicata 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/fornicata
find "$DOWNLOAD_SEQS"/ncbi/fornicata -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/fornicata \;
cat "$DOWNLOAD_SEQS"/ncbi/fornicata/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_fornicata.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/fornicata/

datasets download genome taxon parabasalia 
unzip ncbi_dataset -d "$DOWNLOAD_SEQS"/ncbi/parabasalia
find "$DOWNLOAD_SEQS"/ncbi/parabasalia -type f -name "*.fna" -exec mv {} "$DOWNLOAD_SEQS"/ncbi/parabasalia \;
cat "$DOWNLOAD_SEQS"/ncbi/parabasalia/*.fna > "$DOWNLOAD_SEQS"/ncbi/wgs_parabasalia.fasta
rm -r "$DOWNLOAD_SEQS"/ncbi/parabasalia/
set -e
conda deactivate

#download PR2 database
wget -O "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" "https://github.com/pr2database/pr2database/releases/download/v5.0.0/pr2_version_5.0.0_SSU_taxo_long.fasta.gz"
gzip -d "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta.gz" > "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta"

#filter sequences from PR2 that are not of interest
awk 'BEGIN{IGNORECASE=1} /^>/{printit=1; for (i=1;i<=NF;i++) {if ($i ~ /16S_rRNA|Fungi|uncultured|Arthropoda|Annelida|Craniata|Streptophyta|Chlorophyta/) {printit=0}} } {if (printit) print}' "$DOWNLOAD_SEQS/PR2/pr2_version_5.0.0_SSU_taxo_long.fasta" > "$DOWNLOAD_SEQS/PR2/pr2_clean.fasta"

