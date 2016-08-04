#!/bin/bash

# Damon Polioudakis
# 2016-02-26
# QSUB all samples in input dir to 4_HTSC_Expression.sh which runs HTSeq Counts
# on single-cell RNAseq from VZ and CP demultiplexed and aligned with RNAstarr

# To submit this script:
# qsub -cwd -o logs/4_HTSC_Expression_QSUB_$(date +%Y%m%d).log -e logs/4_HTSC_Expression_QSUB_$(date +%Y%m%d).error -S /bin/bash -V -N HTSC_QSUB_4 -q geschwind.q -l h_data=4G,h_rt=3:00:00 4_HTSC_Expression_QSUB.sh
################################################################################
echo ""
echo "Starting HTSC_Expression_QSUB.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# bam directory contains directories named by sample ID
inParentDir=../data/bam/merged
outDir=../data/htseq/merged
mkdir -p $outDir
################################################################################

for inDir in $inParentDir/HT*; do

  echo ""
  echo "In bam dir:"
  echo ${inDir}
  # sampleID is the directory named by sample ID (removing
  # "../SxaQSEQsXbp060L2_bam/") from the path
  sampleID=$(basename ${inDir})
	echo "In bam sample ID:"
  echo ${sampleID}
  inBam=${inDir}/markdup_sorted.bam
  samDir=${inDir}
  echo "In bam path:"
  echo ${inBam}
  echo "In sam path:"
  echo ${samDir}
  echo "Out dir:"
  echo ${outDir}

  if [ ! -f ${outDir}/${sampleID}_* ]; then

    echo "Submitting ${sampleID} to 4_HTSC_Expression.sh"
    # QSUB 4_HTSC_Expression.sh
    # Note name of job (-N) will only display 1st 10 characters
    outLog=logs/4_HTSC_Expression_$(date +%Y%m%d)_${sampleID}.log
  	outErr=logs/4_HTSC_Expression_$(date +%Y%m%d)_${sampleID}.error

    qsub -cwd -o ${outLog} -e ${outErr} -V -N HT_${sampleID} -S /bin/bash -q geschwind.q -l h_data=36G,h_rt=24:00:00 4_HTSC_Expression.sh ${inBam} ${sampleID} ${outDir} ${samDir}

  else
    echo "HTSC output already exists for:"
    echo ${sampleID}
  fi

done
################################################################################

echo ""
echo "End of 4_HTSC_Expression_QSUB.sh... "$(date)
