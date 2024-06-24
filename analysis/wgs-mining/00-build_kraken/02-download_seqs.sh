#!/bin/sh

# download sequences from wormbase and eupath for building reference genome

#SBATCH --job-name=download_seqs
#SBATCH --nodes=1
#SBATCH --partition medium
#SBATCH --mem=80G
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

set -eu
if [ -f "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env" ]; then
   . "/data/zool-zost/sann7416/island-biogeography-wgs-mining/.env"
fi
# .env DOWNLOAD_SEQS

#wormbase
wget -P "$DOWNLOAD_SEQS"/wormbase  -r -np -A "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

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
