#!/bin/sh

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

# download sequences from wormbase, eupath, ncbi (wgs and nt) and PR2 for building reference databases
set -eu
source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
# .env DOWNLOAD_SEQS, NCBI_DATASETS_CONDA

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
