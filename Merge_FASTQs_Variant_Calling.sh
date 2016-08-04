#!/bin/bash

# Damon Polioudakis
# 2016-04-14
# Select fastqs for human only capture sites and merge lane 1 and 2 to prepare
# for variant calling
# Suggested script call:
#   sh Merge_FASTQs_Variant_Calling.sh 2>&1 | tee logs/Merge_FASTQs_Variant_Calling_$(date +%Y%m%d).log
################################################################################

echo ""
echo "Starting Merge_FASTQs_Variant_Calling.sh"
echo ""
################################################################################

# Define Input Variables and Functions

inSampleID=../analysis/tables/Human_Only_Capture_Sites_10^5Hs_10^5Mm.txt
outDir=../data/fastq/Merged_For_Variant_Calling
mkdir -p ${outDir}
################################################################################

while read sampleID; do
  inFastqsR1=$(ls ../data/fastq/Sxa*/${sampleID}_R1.fastq)
  inFastqsR2=$(ls ../data/fastq/Sxa*/${sampleID}_R2.fastq)
  cat ${inFastqsR1} > ${outDir}/${sampleID}_R1.fastq
  cat ${inFastqsR2} > ${outDir}/${sampleID}_R2.fastq
done < ${inSampleID}
################################################################################

echo ""
echo "End of Merge_FASTQs_Variant_Calling.sh... "$(date)
