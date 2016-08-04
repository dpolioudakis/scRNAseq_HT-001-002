#!/bin/bash

# This script creates a new star index with the new slice junction file
# Use this qsub:
# qsub -cwd -V -N SJindex -S /bin/bash -q geschwind.q -l h_data=16G,h_rt=24:00:00 -pe shared 4  -o logs/Variant_GATK_3_new_index_$(date +%Y%m%d).log -e logs/Variant_GATK_3_new_index_$(date +%Y%m%d).error Variant_GATK_3_new_index.sh

STARcall=/share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR
genomeDir=/geschwindlabshares/RNAseq_singlecellfetal/HT-001-002_DP/source/hg19_2pass ## new directory where you are going to create the genome index
hg19=/geschwindlabshares/RNAseq_singlecellfetal/source/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC.fa
SJ=/geschwindlabshares/RNAseq_singlecellfetal/HT-001-002_DP/data/bam/NoERCC/variant/SJ.all.out.tab  ##this is what is created from 1_RunSJ.sh and 2_SJ.out.ALL.R

## Make new genome file from new splice junction file
mkdir -p $genomeDir
cd $genomeDir

##filter out non-canonical junctions
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' ${SJ} > tmp

##Create Index
$STARcall --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $hg19 --sjdbFileChrStartEnd SJ.all.out.tab.sjdb --sjdbOverhang 49 --runThreadN 4
