#!/bin/bash

# Damon Polioudakis
# 2016-02-24
# Code to submit all samples for FASTQC and STAR Alignment for Stage 1 of
# RNASeq Pipeline

# To submit this script:
# qsub -cwd -o logs/1_FQC_StarAlign_QSUB_NoERCC_$(date +%Y%m%d).log -e logs/1_FQC_StarAlign_QSUB_NoERCC_$(date +%Y%m%d).error -S /bin/bash -V -N SA_QSUB_1 -q geschwind.q -l h_data=4G,h_rt=3:00:00 1_FQC_StarAlign_QSUB_NoERCC.sh

################################################################################
echo ""
echo "Starting 1_FQC_StarAlign_QSUB_NoERCC.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions


# Path to directory that contains fastqs
inParentDir=../data/fastq
# Path for parent directory to output alignment, will create subdirectories
# below for each sample ID
outParentDir=../data/bam/NoERCC
################################################################################

mkdir -p ${outParentDir}
mkdir -p logs

# Alter Sxa* accordingly to match sequencing lane directories
for inDir in ${inParentDir}/Sxa*; do

  # Extract sequencing lane from file path
  seqLane=${inDir##${inParentDir}"/"}
  echo "Sequencing Lane:"
  echo ${seqLane}

  # alter 'R3*_R1_001.fastq.gz' accordingly to match FASTQ files, but make sure that *_R1 is part of wildcard term

  wildC=*R1*fastq
  nFiles=$(ls ${inDir}/${wildC} | wc -l)
  echo "Total files: ${nFiles}"

  for inFastq1 in ${inDir}/${wildC}; do

    echo ""
    echo "fastq Read 1 path:"
    echo ${inFastq1}
    # fastq Read 2 paths
    inFastq2=${inFastq1//R1/R2}
    echo "fastq Read 2 path:"
    echo ${inFastq2}

    # Extract sample ID from file path
    sampleID=${inFastq1##${inDir}"/"}
    sampleID=${sampleID%%[._]R[12]*}
    echo "Sample ID:"
    echo ${sampleID}

    # Designate an output directory for the alignment specific to that sample ID
    outDir=${outParentDir}/${seqLane}/${sampleID}
    echo "Alignment output directory:"
    echo ${outDir}

    if [ ! -d ${outDir} ]; then

      echo "Submitting Star Alignment script for:"
      echo ${inFastq1}
      echo ${inFastq2}

      # File names for QSUB stdout and error output
      outSampLog=logs/1_FQC_StarAlign_NoERCC_$(date +%Y%m%d)_${sampleID}.log
      outSampErr=logs/1_FQC_StarAlign_NoERCC_$(date +%Y%m%d)_${sampleID}.error

      # Note name of job (-N) will only display 1st 10 characters
      qsub -cwd -o ${outSampLog} -e ${outSampErr} -V -S /bin/bash -N A1_${sampleID} -q geschwind.q -l h_data=64G,h_rt=12:00:00 -pe shared 8 1_FQC_StarAlign_NoERCC.sh $outDir $inFastq1 $inFastq2 $inDir $sampleID

    else
      echo "${outDir} already exists..."
    fi

  done

done
################################################################################

echo ""
echo "End of 1_FQC_StarAlign_QSUB_NoERCC.sh... "$(date)
