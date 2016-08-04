#!/bin/bash

# Damon Polioudakis
# 2016-02-24
# Code to format fastq files for Fluidigm demultiplexing
# 0.4_mRNASeqHT_demultiplex.pl script and submit all samples to script

# To submit this script: qsub -cwd -o logs/0.4_mRNASeqHT_demultiplex_QSUB_$(date +%Y%m%d).log -e logs/0.4_mRNASeqHT_demultiplex_QSUB_$(date +%Y%m%d).error -S /bin/bash -V -N FDM_QSUB -q geschwind.q -l h_data=4G,h_rt=3:00:00 0.4_mRNASeqHT_demultiplex_QSUB.sh

# Reminder: make /logs directory in code directory
################################################################################

inParentDir=../data/fastq
# inParentDir=../test/fastq

# Alter Sxa* accordingly to match sequencing lane directories
for inLaneDir in ${inParentDir}/SxaQSEQsXap096L*; do

  echo ""
  echo "In Lane dir:"
  echo ${inLaneDir}

  bnLaneDir=$(basename ${inLaneDir})
  echo "Lane:"
  echo ${bnLaneDir}

  inDir=${inLaneDir}/fastq_not_demultiplexed
  echo "In dir:"
  echo ${inDir}

  echo "Out dir:"
  outDir=${inLaneDir}
  echo ${outDir}

  # Format fastq file names for 0.4_mRNASeqHT_demultiplex.pl
  #   Change . to _ outside of .fastq.gz suffix
  #   Prepend HT_
  #   0.4_mRNASeqHT_demultiplex.pl requires
  #     [ColNum]_[SampleName]_R[1|2].fastq.gz
  # for inFile in ${inDir}/*; do
  #   rnFile=$(echo $(basename $inFile) | awk '/\./{sub(/\./, "_")} {print "HT_"$0}')
  #   echo "Renaming:"
  #   echo "${inFile} to ${inDir}/${rnFile}"
  #   mv ${inFile} ${inDir}/${rnFile}
  # done

  echo "Submitting Fluidigm demultiplex script for:"
  echo ${inDir}

  # File names for QSUB stdout and error output
  outSampLog=logs/0.4_mRNASeqHT_demultiplex_$(date +%Y%m%d)_${bnLaneDir}.log
  outSampErr=logs/0.4_mRNASeqHT_demultiplex_$(date +%Y%m%d)_${bnLaneDir}.error

  # Note name of job (-N) will only display 1st 10 characters
  qsub -cwd -o ${outSampLog} -e ${outSampErr} -V -S /usr/bin/perl -N FDM -q geschwind.q -l h_data=64G,h_rt=12:00:00 -pe shared 8 0.4_mRNASeqHT_demultiplex.pl -i ${inDir}/ -o ${outDir}/

done
