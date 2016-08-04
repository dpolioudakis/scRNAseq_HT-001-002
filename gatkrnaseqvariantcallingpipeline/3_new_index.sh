## Use this qsub: qsub -cwd -V -N SJindex -S /bin/bash -q geschwind.q -l h_data=16G,h_rt=24:00:00 -pe shared 4 3_new_index.sh
#!/bin/bash

##This script creates a new star index with the new slice junction file

STARcall=/share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR
genomeDir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/hg19_2pass ## new directory where you are going to create the genome index
hg19=/geschwindlabshares/eQTL/GenomeBuild/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa 
SJ=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/hg19_2pass_R2/SJ.all.out.tab  ##this is what is created from 1_RunSJ.sh and 2_SJ.out.ALL.R 



## Make new genome file from new splice junction file
cd $genomeDir

##filter out non-canonical junctions
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' SJ.all.out.tab > SJ.all.out.tab.sjdb

##Create Index
$STARcall --runMode genomeGenerate --genomeDir $genomeDir1 --genomeFastaFiles $hg19 --sjdbFileChrStartEnd $SJ.sjdb --sjdbOverhang 75 --runThreadN 4
