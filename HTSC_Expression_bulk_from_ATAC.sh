#!/bin/bash

# Damon Polioudakis
# 2016-02-08
# Run HTSeq Counts on bulk RNAseq VZ and CP from Luis and Jason ATAC samples
################################################################################
echo "Starting HTSC_Expression_bulk_from_ATAC.sh... "$(date)
################################################################################

# Define Input Variables and Functions

inBam=$1
sampleID=$2
outDir=$3
samDir=$4

gtf=../source/gencode.v19.annotation.gtf
outHTSC=${outDir}/${sampleID}_exon_union_count.txt
# Samtools sort output path, samtools adds .bam to ouput so leave suffix off # file name
bamNameSrtdPfx=${outDir}/tmp_${sampleID}_name_sorted
# Samtools sort output path
bamNameSrtd=${outDir}/tmp_${sampleID}_name_sorted.bam
################################################################################

echo "Started HTSC exon on...: "${inBam}
# Unsure why samtools needs an output file name (tmp.bam) in order to output to standard out, no output file is actually written, but the standard out can then be piped to samtools view
if [ ! -s outHTSC ]; then

  mkdir -p $outDir

  echo "In bam:"
  echo ${inBam}

  # Samtools sort by bam by name (required by HTSeq)
  echo "Start samtools sorting by name on...:"
  echo ${inBam}
  # samtools adds .bam to ouput so left off .bam suffix from variable
  samtools sort -n ${inBam} ${bamNameSrtdPfx}
  echo "Done samtools sorting by name, output:"
  echo ${bamNameSrtd}

  # Convert name sorted bam to sam and pipe to HTseq
  echo "Started HTSC exon on...:"
  echo ${bamNameSrtd}
  # "-" tells HTseq to input from pipe
  samtools view -h ${bamNameSrtd} | /share/apps/anaconda/bin/htseq-count --stranded=no --mode=union --type=exon - $gtf >> ${outHTSC}

  rm -f ${bamNameSrtd}
  echo "Saving HTSC exon results to...: "${outHTSC}

else
  echo "${outHTSC} already exists"
fi

date +%F%t%T
echo "End of HTSC_Expression.sh..."
