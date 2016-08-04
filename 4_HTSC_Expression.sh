#!/bin/bash

# Damon Polioudakis
# 2016-02-26
# Run HTSeq Counts
################################################################################
echo ""
echo "Starting 4_HTSC_Expression.sh... "$(date)
################################################################################

inBam=$1
sampleID=$2
outDir=$3
samDir=$4

gtf=../../source/gencode.v19.annotation.gtf
outHTSCeu=${outDir}/${sampleID}_exon_union_count.txt
outHTSCgu=${outDir}/${sampleID}_gene_union_count.txt
# Samtools sort output path, samtools adds .bam to ouput so leave suffix off
# file name
bamNameSrtdPfx=${samDir}/markdup_name_sorted
# Samtools sort output path
bamNameSrtd=${samDir}/markdup_name_sorted.bam
# If Samtools view conversion of .bam to .sam pipe into htseq-count does not
# work, can write .sam to disk and then input into htseq-count
# samPath=${samDir}/PEmatched_markdup_name_sorted.sam
################################################################################

echo ""
echo "In bam:"
echo ${inBam}
echo "SampleID:"
echo ${sampleID}
echo "Out dir:"
echo ${outDir}

# Samtools sort by name
if [ ! -s ${bamNameSrtd} ]; then
  # Samtools sort by bam by name (required by HTSeq)
  echo "Start samtools sorting by name on...:"
  echo ${inBam}
  # samtools adds .bam to ouput so left off .bam suffix from variable
  samtools sort -n ${inBam} ${bamNameSrtdPfx}
  echo "Done samtools sorting by name, output:"
  echo ${bamNameSrtd}
else
  echo "${bamNameSrtd} already exists"
fi

# HSTC exon union
if [ ! -s ${outHTSCeu} ]; then
  # Convert name sorted bam to sam and pipe to HTseq
  echo "Started HTSC on...:"
  echo ${bamNameSrtd}
  # "-" tells HTseq to input from pipe
  samtools view -h ${bamNameSrtd} | /share/apps/anaconda/bin/htseq-count --stranded=no --mode=union --type=exon - $gtf >> ${outHTSCeu}
  echo "Saving HTSC exon union results to...: "${outHTSCeu}
else
  echo "${outHTSCeu} already exists"
fi

# HSTC gene union
if [ ! -s ${outHTSCgu} ]; then
  # HSTC exon union
  # "-" tells HTseq to input from pipe
  samtools view -h ${bamNameSrtd} | /share/apps/anaconda/bin/htseq-count --stranded=no --mode=union --type=gene - $gtf >> ${outHTSCgu}

  # Remove name sorted bam after HTseq Counts runs both exon and gene union
  rm -f ${bamNameSrtd}
  echo "Saving HTSC gene union results to...: "${outHTSCgu}
else
  echo "${outHTSCgu} already exists"
fi
################################################################################

echo ""
echo "End of 4_HTSC_Expression.sh... "$(date)
