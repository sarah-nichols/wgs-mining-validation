#!/bin/sh

#SBATCH --job-name=filter_host
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --time=130:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk
#SBATCH --array=1-72


module purge
module load SAMtools/1.14-GCC-11.2.0
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

source "/data/biol-bird-parasites/sann7416/wgs-mining-validation/src/.env"
#.env containing SAMPLE_IDS_ZOST, BAM_ZOST_PATH, HOST_ZOST_HSS, HOST_ZOST_IMG, FILTER_METRICS_ZOST 

# Get the $SLURM_ARRAY_TASK_ID-th line from the file
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_IDS_ZOST")

BAM_OUT="$BAM_ZOST_PATH/${BAM##/*/}"
echo "Copying to $BAM_OUT"
cp "$BAM" "$BAM_OUT"
echo "Editing header with samtools"
BAM_HEADER="$BAM_ZOST_PATH/header_${BAM##/*/}"
samtools view -H "$BAM_OUT" | sed 's/,//' | samtools reheader - "$BAM_OUT" > "$BAM_HEADER"
rm "$BAM_OUT"
BAM_FILTERED_PAIRED="$BAM_ZOST_PATH/paired/${BAM##/*/}"
echo "Paired output: $BAM_FILTERED_PAIRED"
BAM_FILTERED_UNPAIRED="$BAM_ZOST_PATH/unpaired/${BAM##/*/}"
echo "Unpaired output: $BAM_FILTERED_UNPAIRED"
gatk --java-options "-Xmx80G" PathSeqFilterSpark  \
--input "$BAM_HEADER" \
--paired-output "$BAM_FILTERED_PAIRED" \
--unpaired-output "$BAM_FILTERED_UNPAIRED" \
--min-adapter-length 1 \
--min-clipped-read-length 60 \
--kmer-file "$HOST_ZOST_HSS" \
--is-host-aligned TRUE \
--filter-bwa-image "$HOST_ZOST_IMG" \
--filter-metrics "$FILTER_METRICS_ZOST/metrics_${BAM##/*/}" \
--bam-partition-size 4000000
rm "$BAM_HEADER"


 
