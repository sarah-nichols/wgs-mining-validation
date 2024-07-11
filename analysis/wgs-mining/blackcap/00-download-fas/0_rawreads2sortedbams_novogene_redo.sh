#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --array=1-2:1
#SBATCH --job-name=Filter2Bam_Pipeline
#SBATCH --output=Filter2Bam.%A_%a.out
#SBATCH --error=Filter2Bam.%A_%a.error

#######################################################################################################
# SET SAMPLE NAME AND DIRECTORY
# This pipeline requires the number of jobs to match the number of lines
# in file sample.list.txt (update #SBATCH --array=1-137:1)
#######################################################################################################

# STEP 1:
# We will need to create a text file specifying the name of samples
# we want to process and the directory that raw reads are stored in,
# using 3 samples as an exaple, this will look like this:

# NOR140 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR141 /data/zool-zost/Novogene/NorfolkIsland/NOR142
# NOR142 /data/zool-zost/Novogene/NorfolkIsland/NOR142

Sample_List=/data/zool-zost/sjoh4959/projects/0.0_gwas_migration/data/resources/raw_reads_novogene4

# STEP 2:
# Use slurm array task ID to alocate sample name and directory
SAMPLE_NAME=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $1}')
SAMPLE_DIRECTORY=$(cat $Sample_List | head -n $SLURM_ARRAY_TASK_ID | tail -1 | awk {'print $2}')

# STEP 3:
#move into sample directory
cd $SAMPLE_DIRECTORY

#######################################################################################################
# CONDUCT FILTERING OF WGS RAW READS
# We will use FastP - a tool designed to provide fast all-in-one preprocessing for FastQ files
# https://github.com/OpenGene/fastp
#######################################################################################################

# STEP 1:
# Define path to Fastp:
FASTP=/data/zool-zost/BIN/fastp

# STEP 2:
# Set up for loop to conduct filtering for each read pair
for ReadPair in `ls ${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4 -d'_'`
do
  
  #Use Fastp to conduct automated filtering of fastq files
  #Note: based on initial test we will trim the first 10bp from start of each read
  $FASTP \
  -i ${ReadPair}_1.fq.gz \
  -o Filtered_${ReadPair}_1.fq.gz \
  -I ${ReadPair}_2.fq.gz \
  -O Filtered_${ReadPair}_2.fq.gz \
  --trim_front1 10 \
  --trim_front2 10

   #Rename QC report and move to fastp qc report folder
  mv fastp.html /data/zool-zost/sjoh4959/projects/0.0_gwas_migration/reports/raw_read_QC/${ReadPair}.html
  mv fastp.json /data/zool-zost/sjoh4959/projects/0.0_gwas_migration/reports/raw_read_QC/${ReadPair}.json

done


#######################################################################################################
# MAP FILTERED READS TO REFERENCE ASSEMBLY
#######################################################################################################

# STEP 1:
# Define path to BWA and reference assembly
# Note: If not already done, will need to index reference assembly ($BWA index $REF)

BWA=/data/zool-zost/BIN/bwa/bwa
REF=/data/zool-zost/sjoh4959/projects/0.0_gwas_migration/data/Zlat_2_Tgutt_ZW.fasta

# STEP 2:
# Load required modules in ARC
ml SAMtools/1.14-GCC-11.2.0

# STEP 3:
# For each pair of reads conduct mapping using BWA MEM
for ReadPair in `ls Filtered_${SAMPLE_NAME}_*_1.fq.gz | cut -f1,2,3,4,5 -d'_'`
do
  	$BWA mem $REF \
  	${ReadPair}_1.fq.gz \
  	${ReadPair}_2.fq.gz	\
  	-R "@RG\tID:${ReadPair}\tSM:${SAMPLE_NAME}" \
  	| samtools view -bS - \
  	| samtools sort - > ${ReadPair}.bam
done


#######################################################################################################
# MERGE SAMPLE BAMS INTO A SINGLE FILE
# Currently we have a bam file for each read pair per sample. We will now merge these to create a
# single bam file per sample
#######################################################################################################

# Count number of bam files per sample
BAM_Count=$(ls Filtered_${SAMPLE_NAME}*.bam | wc -l)

# If number of bams is greater than 1 then merge bams into a single file using samtools merge
# else, rename single bam file
if [ $BAM_Count -gt 1 ]
then
  samtools merge ${SAMPLE_NAME}.bam Filtered_${SAMPLE_NAME}*.bam -f
# rm Filtered_${SAMPLE_NAME}*.bam
else
  mv Filtered_${SAMPLE_NAME}*.bam ${SAMPLE_NAME}.bam
fi

samtools sort ${SAMPLE_NAME}.bam > ${SAMPLE_NAME}.sorted.bam
#rm ${SAMPLE_NAME}.bam

samtools index ${SAMPLE_NAME}.sorted.bam

#######################################################################################################
# REMOVE FILTERED READ FILES
# As we have storage issues we will remove filtered read files to free up space
#######################################################################################################

#rm Filtered_${SAMPLE_NAME}_*.fq.gz
rm *.log
