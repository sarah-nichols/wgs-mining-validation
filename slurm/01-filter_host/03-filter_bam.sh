#!/bin/sh

#SBATCH --job-name=filter_host
#SBATCH --nodes=1
#SBATCH --partition long
#SBATCH --mem=80G
#SBATCH --clusters=all
#SBATCH --time=130:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user sarah.nichols@biology.ox.ac.uk
#SBATCH --array=1-117


module purge
module load SAMtools/1.14-GCC-11.2.0
module load GATK/4.1.5.0-GCCcore-9.3.0-Java-1.8

set -x
if [ -f "/data/zool-zost/sann7416/wgs-mining-validation/.env"]; then
   .  "/data/zool-zost/sann7416/wgs-mining-validation/.env"
fi

# Get the $SLURM_ARRAY_TASK_ID-th line from the file
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_IDS")

BAM_OUT="$OUTPATH/${BAM##/*/}"
echo "Copying to $BAM_OUT"
cp "$BAM" "$BAM_OUT"
echo "Editing header with samtools"
BAM_HEADER="$OUTPATH/header_${BAM##/*/}"
samtools view -H "$BAM_OUT" | sed 's/,//' | samtools reheader - "$BAM_OUT" > "$BAM_HEADER"
rm "$BAM_OUT"
BAM_FILTERED_PAIRED="$OUTPATH/paired_${BAM##/*/}"
BAM_FILTERED_UNPAIRED="$OUTPATH/unpaired_${BAM##/*/}"
gatk --java-options "-Xmx80G" PathSeqFilterSpark  \
--input "$BAM_HEADER" \
--paired-output "$BAM_FILTERED_PAIRED" \
--unpaired-output "$BAM_FILTERED_UNPAIRED" \
--min-adapter-length 1 \
--min-clipped-read-length 60 \
--kmer-file "$HOST_HSS" \
--is-host-aligned TRUE \
--filter-bwa-image "$HOST_IMG" \
--filter-metrics "$FILTER_METRICS/metrics_${BAM##/*/}" \
--bam-partition-size 4000000
rm "$BAM_HEADER"


 
