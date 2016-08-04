#!/bin/bash

# Damon Polioudakis
# 2016-02-08
# QSUB all samples listed in ../metadata/VZCP_sampleinfo.csv BAM column to HTSC_Expression_bulk_from_ATAC.sh which runs HTSeq Counts
# Samples are bulk RNAseq from Luis and Jason ATAC samples
# To submit this script:
# qsub -cwd -o logs/HTSC_Expression_bulk_from_ATAC_QSUB_$(date +%Y%m%d).log -e logs/HTSC_Expression_bulk_from_ATAC_QSUB_$(date +%Y%m%d).error -S /bin/bash -V -N HTSC_QSUB_4 -q geschwind.q -l h_data=4G,h_rt=3:00:00 HTSC_Expression_bulk_from_ATAC_QSUB.sh

echo "Starting HTSC_Expression_bulk_from_ATAC_QSUB.sh..."

outDir=../data/htseq/bulk_VZ_CP_from_ATAC
mkdir -p $outDir

for inBam in $(awk 'FS="," {print $20}' < ../metadata/VZCP_sampleinfo.csv); do

  echo "In bam: "${inBam}
  # sampleID is the directory named by sample ID (removing
  # "../SxaQSEQsXbp060L2_bam/") from the path
  sampleID=$(basename ${inBam} .merge.remarkdup.bam)
	echo "In bam sample name: "${sampleID}

  # QSUB 4_HTSC_Expression.sh
  # Note name of job (-N) will only display 1st 10 characters
  outLog=logs/HTSC_Expression_bulk_from_ATAC_$(date +%Y%m%d)_${sampleID}.log
  outErr=logs/HTSC_Expression_bulk_from_ATAC_$(date +%Y%m%d)_${sampleID}.error

  qsub -cwd -o ${outLog} -e ${outErr} -V -N HT_${sampleID} -S /bin/bash -q geschwind.q -l h_data=36G,h_rt=24:00:00 HTSC_Expression_bulk_from_ATAC.sh ${inBam} ${sampleID} ${outDir}

done

date +%F%t%T
echo "End of HTSC_Expression_QSUB.sh..."
