#!/bin/bash

# Damon Polioudakis
# 2016-03-01
# Count reads found in both bam files

# Suggest script call:
# sh Count_Reads_In_Both_Bams.sh 2>&1 | tee logs/Count_Reads_In_Both_Bams_$(date +%Y%m%d).log
# Reminder: make /logs directory in code directory

# Outputs counts from all lanes into 1 file
################################################################################
echo ""
echo "Starting Count_Reads_In_Both_Bams.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# Path to directory that contains bams
inParentDir=../data/bam
# Path for parent directory to output processing files, will create
# subdirectories for each sample ID
outParentDir=../data/count_reads_in_human_and_mouse
# Output compiled table of numbers of reads aligning to both human and mouse for
# each sample
outStats=${outParentDir}/Reads_Aligned_Human_And_Mouse_Stats.txt
################################################################################

mkdir -p ${outParentDir}

# Header for output compiled statistics file
echo -e "SampleID\tReads_In_Human\tReads_In_Mouse\tAlign_To_Both\tAlign_To_Human_Only\tAlign_To_Mouse_Only" > ${outStats}

# Convert bams to sams, extract read names, and sort unique
for inBamDir in ${inParentDir}/*/SxaQSEQsXap096L*/*; do

  echo ""
  echo "Converting bam to sam, extracting read names, and sort unique..."

  inBam=${inBamDir}/Aligned.sortedByCoord.out.bam
  echo "In Bam:"
  echo ${inBam}

  outDir=${outParentDir}/${inBamDir##${inParentDir}/}
  echo "Out Dir:"
  echo ${outDir}
  mkdir -p ${outDir}

  # -q 254 option selects reads with map quality score > 254
  # The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads
  samtools view -q 254 ${inBam} | cut -f1 | sort -u > ${outDir}/Reads.txt
done

# Compare sorted read name files line by line to identify reads common to both
# mouse and human alignments
for outDirHuman in ${outParentDir}/NoERCC/Sxa*/*; do

  echo ""
  echo "Identifying reads common to both bams..."

  echo "Out Dir Human:"
  echo ${outDirHuman}

  outDirMouse=$(echo ${outDirHuman} | sed s/NoERCC/mouse/)
  echo "Out Dir Mouse:"
  echo ${outDirMouse}

  outDir=${outParentDir}/${outDirHuman##../data/compare_human_to_mouse/NoERCC/}
  echo "Out Dir:"
  echo ${outDir}
  mkdir -p ${outDir}

  sampleID=$(basename ${outDirHuman})
  echo "Sample ID:"
  echo ${sampleID}

  # Reads Aligned to Both
  comm -12 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Both.txt

  # Reads Aligned to Human Only
  comm -23 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Human_Only.txt

  # Reads Aligned to Mouse Only
  comm -13 ${outDirHuman}/Reads.txt ${outDirMouse}/Reads.txt > ${outDir}/Align_To_Mouse_Only.txt

  # Compile and write out numbers aligning to both, human only, or mouse only
  echo -ne ${sampleID}"\t" >> ${outStats}
  echo -ne $(wc -l ${outDirHuman}/Reads.txt | cut -f1 -d" ")"\t" >> ${outStats}
  echo -ne $(wc -l ${outDirMouse}/Reads.txt | cut -f1 -d" ")"\t" >> ${outStats}
  echo -ne $(wc -l ${outDir}/Align_To_Both.txt | cut -f1 -d" ")"\t" >> ${outStats}
  echo -ne $(wc -l ${outDir}/Align_To_Human_Only.txt | cut -f1 -d" ")"\t" >> ${outStats}
  echo -e $(wc -l ${outDir}/Align_To_Mouse_Only.txt | cut -f1 -d" ") >> ${outStats}

done

# Remove intermediate files
rm -r ${outParentDir}/mouse
rm -r ${outParentDir}/NoERCC
rm -r ${outParentDir}/SxaQSEQsXap096L*
################################################################################

echo ""
echo "End of Count_Reads_In_Both_Bams.sh... "$(date)
