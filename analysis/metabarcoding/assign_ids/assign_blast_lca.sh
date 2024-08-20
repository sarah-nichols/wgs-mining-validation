#!/bin/bash

#SBATCH --job-name=blast2lca
#SBATCH --output=blast2lca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -o blast2lca_%a.out # Standard output
#SBATCH -e blast2lca_%a.err # Standard error
#SBATCH --partition medium
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=6:00:00
#SBATCH --clusters=ALL
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk

source /data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env

module load Anaconda3/2024.02-1


USAGE="Usage: $(basename "$0") -B <blast percent ID value 0-100> -T <blast2lca top percent 1-10> -D <abs path of megan nucl database> \n
The script takes the output from blast (located in the blast_out directory), and applies the following: \n
step1: Merge results if blast was run in array mode (automatically detected);
step2: Filter blast results by percentage ID (0-100). The user provides a 
       minimum percentage ID  with -B for filtering blast results. We recommend 85-95;
step3: Run the megan2lca (lowest common ancestor) algorithm to determine taxonomic likelihood of 
       ASV at all taxonomic levels; -T should be set to a value between 1 and 10
step4: Formatting of the output of blast2lca for downstream analyses
       This includes taxon path files and a summary file
Here is an example of how to run the script: \n
qsub b2m_scripts/02_run_blast2lca.sh -B 90 -T 2 -D /shared/genomicsdb2/shared/megan/megan-nucl-Feb2022.db \n
The script assumes blast results (whether from simple or array mode) are located in the directory 
blast_out, and also saves intermediate and final files to the same directory. \n\n"

## List arguments
while getopts B:T:D: flag; do
	case "${flag}" in
		B) BPI=${OPTARG};;
		T) TOP_PER=${OPTARG};;
		D) MEG_DB=${OPTARG};;
	esac
done

## Check mandatory arguments
shift $((OPTIND-1))
if [ -z "${BPI}" ] || [ -z "${TOP_PER}" ] || [ -z "${MEG_DB}" ]; then
   printf "\n\n${USAGE}" >&2; exit 1
fi

echo "blast percent identity value = " ${BPI}
echo "blast2lca top percent value = "${TOP_PER}


## Define path variables
MAIN_DIR=$PWD
OUT_DIR="wgs-mining-validation/data/metabarcoding/blast"

## Step 1 : Check if run in array mode and, if so, merge the chunks to create all_blast.out.tab
if ls ${MAIN_DIR}/${OUT_DIR}/blast_results_*.out 1> /dev/null 2>&1; then
  echo "Blast was run in array mode, merging chunks..."
  cat ${MAIN_DIR}/${OUT_DIR}/blast_results_*.out > ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab
fi


## Step 2: Remove additional taxa information and filter by user specified blast percentage identity (BPI)
cut -f1-12 ${MAIN_DIR}/${OUT_DIR}/all_blast.out.tab | awk -v var="${BPI}" '$3 >= var' > ${MAIN_DIR}/${OUT_DIR}/filtered_blast.out.tab

source activate $MEGAN_ENV

export LCA_PATH=/data/biol-bird-parasites/sann7416/conda_environments/megan_env/bin/blast2lca
export PATH=/data/biol-bird-parasites/sann7416/conda_environments/megan_env/bin/:$PATH

## Step 3: Run the Megan blast2lca (lowest common ancestor) algorithm with a user specified top percent threshold (TOP_PER)
/data/biol-bird-parasites/sann7416/conda_environments/megan_env/bin/blast2lca -i ${MAIN_DIR}/${OUT_DIR}/filtered_blast.out.tab -m BlastN -o ${MAIN_DIR}/${OUT_DIR}/megan_full_out.tsv -mdb ${MEG_DB} -top ${TOP_PER}

## Step 4: Take the full megan output and produce summary files for downstream analysis.
## This is actually 5 separate "one-liners": piped awk (and sed) commands to format the megan_full_out.tsv file.
##
## The first two commands produce temporary files (tmp_step*) that will be deleted at the end of the script. These are necessary
## because the number of fields (columns) can vary in the megan output.
##
## The third command produces a taxon path file containing the numbers of blast hits assigned at each taxonomic level.
## In most cases this won't be used, but can be diagnisotic if many ASVs have only higher taxonomic assignment.
##
## The fourth command produces a taxon path file where assignment is provided up to the point at which 100% of blast
## hits match (the lowest common ancestor). Lower taxonomic levels are assigned values of NA.
##
## The fifth command produces a taxon summary file with just the lowest common ancestor and it's taxonomic rank.

## Step 4 Command 1
sed -e 's/;/,/g' ${MAIN_DIR}/${OUT_DIR}/megan_full_out.tsv | sed -e 's/s__/b__/3' | \
    awk -F ',' '{if ($0 ~ /d__/) {print $0} \
                 else {print $0 ",d__missing, NA"}}' | \
    awk -F ',' '{if ($0 ~ /k__/) {print $0} \
                 else if ($0 !~ /d__Eukaryota/) {split ($3, DOM, "__"); print $0 ",k__" DOM[2] ", NA"} \
                 else {print $0 ",k__missing, NA"}}'| 
    awk -F ',' '{if ($0 !~ /p__/) {print $0 ",p__missing, NA"} \
                 else if ($0 ~ /p__.+p__/ || $0 !~ /p__unknown/) {print $0} \
                 else {print $0 ",p__missing, NA"}}' | \
    awk -F ',' '{if ($0 !~ /c__/) {print $0 ",c__missing, NA"} \
                 else if ($0 ~ /c__.+c__/ || $0 !~ /c__unknown/) {print $0} \
                 else {print $0 ",c__missing, NA"}}' | \
    awk -F ',' '{if ($0 !~ /o__/) {print $0 ",o__missing, NA"} \
                 else if ($0 ~ /o__.+o__/ || $0 !~ /o__unknown/) {print $0} \
                 else {print $0 ",o__missing, NA"}}' | \
    awk -F ',' '{if ($0 !~ /f__/) {print $0 ",f__missing, NA"} \
                 else if ($0 ~ /f__.+f__/ || $0 !~ /f__unknown/) {print $0} \
                 else {print $0 ",f__missing, NA"}}' | \
    awk -F ',' '{if ($0 !~ /g__/) {print $0 ",g__missing, NA"} \
                 else if ($0 ~ /g__.+g__/ || $0 !~ /g__unknown/) {print $0} \
                 else {print $0 ",g__missing, NA"}}' | \
    awk -F ',' '{if ($0 !~ /s__/) {print $0 ",s__missing, NA"} \
                 else if ($0 ~ /s__.+s__/ || $0 !~ /s__unknown/) {print $0} \
                 else {print $0 ",s__missing, NA"}}' | \
    awk -F ',' '{if ($0 ~ /b__/ || $0 ~ /v__/) {print $0} \
                 else {print $0 ",b__NA, NA"}}' \
> ${MAIN_DIR}/${OUT_DIR}/tmp_step1_megan_prep.tsv

## Step 4 Command 2
awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /d__/ && $i !~ /d__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' ${MAIN_DIR}/${OUT_DIR}/tmp_step1_megan_prep.tsv | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /k__/ && $i !~ /k__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /p__/ && $i !~ /p__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /c__/ && $i !~ /c__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /o__/ && $i !~ /o__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /f__/ && $i !~ /f__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /g__/ && $i !~ /g__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if ($i ~ /s__/ && $i !~ /s__unknown/) {print $1 ";" $i ";" $(i+1) ",\t" $0} } }' | \
    awk -F ',' '{ for (i=1; i<=NF; ++i) { if (($i ~ /b__/ && $i !~ /b__unknown/) ||($i ~ /v__/ && $i !~ /v__unknown/)) {print $1 ";" $i ";" $(i+1) } } }' \
> ${MAIN_DIR}/${OUT_DIR}/tmp_step2_megan_prep.tsv

## Step 4 Command 3
awk -F ';' '{if ($4 ~ /k__missing/) {split ($2, DOM, "__"); print $1 ";" $2 ";" $3 ";k__unknown " DOM[2] " kingdom; NA; "$6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else {print $0}}' ${MAIN_DIR}/${OUT_DIR}/tmp_step2_megan_prep.tsv | \
    awk -F ';' '{if ($6 ~ /p__missing/ && $4 !~ /k__unknown/) {split ($4, KIN, "__"); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";p__unknown " KIN[2] " phylum; NA;" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else if ($6 ~ /p__missing/ && $4 ~ /k__unknown/) {split ($4, KIN, " "); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";p__unknown " KIN[2] " phylum; NA;" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19}
                 else {print $0}}' | \
    awk -F ';' '{if ($8 ~ /c__missing/ && $6 !~ /p__unknown/) {split ($6, PHY, "__" ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";c__unknown " PHY[2] " class; NA;" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else if ($8 ~ /c__missing/ && $6 ~ /p__unknown/) {split ($6, PHY, " " ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";c__unknown " PHY[2] " class; NA;" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else {print $0}}' | \
    awk -F ';' '{if ($10 ~ /o__missing/ && $8 !~ /c__unknown/) {split ($8, CLA, "__" ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";o__unknown " CLA[2] " order; NA;" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else if ($10 ~ /o__missing/ && $8 ~ /c__unknown/) {split ($8, CLA, " " ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";o__unknown " CLA[2] " order; NA;" $12 ";" $13 ";" $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else {print $0}}' | \
    awk -F ';' '{if ($12 ~ /f__missing/ && $10 !~ /o__unknown/) {split ($10, ORD, "__" ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";f__unknown " ORD[2] " family; NA;"  $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else if ($12 ~ /f__missing/ && $10 ~ /o__unknown/) {split ($10, ORD, " " ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";f__unknown " ORD[2] " family; NA;"  $14 ";" $15 ";" $16 ";" $17 ";" $18 ";" $19} \
                 else {print $0}}' | \
    awk -F ';' '{if ($14 ~ /g__missing/ && $12 !~ /f__unknown/) {split ($12, FAM, "__" ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";g__unknown " FAM[2] " genus; NA;"  $16 ";" $17 ";" $18 ";" $19} \
                 else if ($14 ~ /g__missing/ && $12 ~ /f__unknown/) {split ($12, FAM, " " ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";g__unknown " FAM[2] " genus; NA;"  $16 ";" $17 ";" $18 ";" $19} \
                 else {print $0}}' | \
    awk -F ';' '{if ($16 ~ /s__missing/ && $14 !~ /g__unknown/) {split ($14, GEN, "__" ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";s__unknown " GEN[2] " species; NA;" $18 ";" $19} \
                 else if ($16 ~ /s__missing/ && $14 ~ /g__unknown/) {split ($14, GEN, " " ); print $1 ";" $2 ";" $3 ";" $4 ";" $5 ";" $6 ";" $7 ";" $8 ";" $9 ";" $10 ";" $11 ";" $12 ";" $13 ";" $14 ";" $15 ";s__unknown " GEN[2] " species; NA;" $18 ";" $19} \
                 else {print $0}}' \
> ${MAIN_DIR}/${OUT_DIR}/megan_taxonpath_withcounts.tsv

## Note that we assign an environmental variable, "HITS", which is the number of blast hits that match at a given taxonomic level.
## To confidently and correctly assign taxa this should be set to 100(%), so we do not advise changing the value. However, it can be 
## reduced to increase taxonomic assignment at lower levels. Just be aware that this will likely produce spurious or arbitrary results.

HITS=100

## Step 4 Command 4
awk -v var="${HITS}" -F ';' '
{ if (NF == 19) { \
        if ($(NF) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" $(NF-9) "\t" $(NF-7) "\t" $(NF-5) "\t" $(NF-3) "\t" $(NF-1)} \
        else if ($(NF-2) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" $(NF-9) "\t" $(NF-7) "\t" $(NF-5) "\t" $(NF-3) "\t" "x__NA"} \
        else if ($(NF-4) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" $(NF-9) "\t" $(NF-7) "\t" $(NF-5) "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-6) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" $(NF-9) "\t" $(NF-7) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-8) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" $(NF-9) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-10) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" $(NF-11) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-12) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" $(NF-13) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-14) >= var) {print $1 "\t" $(NF-17) "\t" $(NF-15) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else if ($(NF-16) >= var) {print $1 "\t" $(NF-17) "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
        else {print $1 "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
    } \
    else {print $1 "\t" "z__check_megan_full_out" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA" "\t" "x__NA"} \
}' ${MAIN_DIR}/${OUT_DIR}/megan_taxonpath_withcounts.tsv | sed -e 's/.__//g' | tr ' ' '_' \
> ${MAIN_DIR}/${OUT_DIR}/megan_taxonpath_out.tsv

## Step 4 Command 5
awk -v var="${HITS}" -F ';' '
{ \
    if (NF == 19) { \
        if ($(NF) >= var) {print $1 "\t" $(NF-1)} \
        else if ($(NF-2) >= var) {print $1 "\t" $(NF-3)} \
        else if ($(NF-4) >= var) {print $1 "\t"  $(NF-5)} \
        else if ($(NF-6) >= var) {print $1 "\t"  $(NF-7)} \
        else if ($(NF-8) >= var) {print $1 "\t" $(NF-9)} \
        else if ($(NF-10) >= var) {print $1 "\t" $(NF-11)} \
        else if ($(NF-12) >= var) {print $1 "\t" $(NF-13)} \
        else if ($(NF-14) >= var) {print $1 "\t" $(NF-15)} \
        else if ($(NF-16) >= var) {print $1 "\t" $(NF-17)} \
        else {print $1 "\t" "x__NA"} \
    } \
    else {print $1 "\t" "y__check_megan_full_out"} \
}' ${MAIN_DIR}/${OUT_DIR}/megan_taxonpath_withcounts.tsv | awk -v OFS="\t" -F '[\t__]' '{print $1 "_" $2 "\t" $3 "\t" $5}' | tr ' ' '_' \
> ${MAIN_DIR}/${OUT_DIR}/megan_summary_out.tsv

# Remove tmp files now
 rm ${MAIN_DIR}/${OUT_DIR}/tmp_step*_megan_prep.tsv
 
 conda deactivate