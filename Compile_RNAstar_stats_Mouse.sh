#!/bin/bash

# Damon Polioudakis
# 2015-02-18
# Compile RNAstar stats from RNAstar Log.final.out
# Suggested script call:
#   sh Compile_RNAstar_stats_Mouse.sh 2>&1 | tee logs/Compile_RNAstar_stats_Mouse_$(date +%Y%m%d).log
################################################################################
echo ""
echo "Starting Compile_RNAstar_stats_Mouse.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# bam directory contains directories named by sample ID
inParentDir=../data/bam/mouse
outDir=../metadata
mkdir -p $outDir
################################################################################

for inLane in ${inParentDir}/Sxa*; do
  echo ""
  echo "#####"
  echo ""
  echo "Lane:"
  echo ${inLane}
  bnLane=$(basename ${inLane})
  echo "Basename:"
  echo ${bnLane}
  outFile=${outDir}/RNAstar_Stats_Mouse_${bnLane}.txt

  for inDir in ${inLane}/*; do
    echo ""
    echo "In dir:"
    echo ${inDir}

    # Sample name is directory name
    bnInDir=$(basename ${inDir})
    echo "Basename:"
    echo ${bnInDir}

    # Path to RNA Star log file
    inLog=${inDir}/Log.final.out
    echo "Log:"
    echo ${inLog}

    # If output compiled stats file does not exit, make and add header
    # Else append stats to file
    if [ ! -f ${outFile} ]; then

      # Add SampleID column
      echo -ne "SampleID\t" > ${outFile}
      # Header from 1st log file read
      # Remove all consecutive tabs and spaces from beginning and end of field
      # and remove | from end of field
      awk '{FS="\t"; ORS="\t"} {gsub(/^[ \t]+/,"",$1); gsub(/[ \t\|]+$/,"",$1); if(NR>5) print $1}' ${inLog} >> ${outFile}
      echo "" >> ${outFile}

      # Stats from 1st log file read
      # Add SampleID
      echo -ne ${bnInDir}"\t" >> ${outFile}
      # Stats
      awk '{FS="\t"; ORS="\t"} {if(NR>5) print $2}' ${inLog} >> ${outFile}
      echo "" >> ${outFile}
    else
      # Stats from other log files
      # Add SampleID
      echo -ne ${bnInDir}"\t" >> ${outFile}
      # Stats
      awk '{FS="\t"; ORS="\t"} {if(NR>5) print $2}' ${inLog} >> ${outFile}
      echo "" >> ${outFile}
    fi
  done
done
################################################################################

echo ""
echo "End of Compile_RNAstar_stats_Mouse.sh... "$(date)
