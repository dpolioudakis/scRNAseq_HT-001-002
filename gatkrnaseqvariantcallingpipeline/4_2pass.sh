## Use this qsub: qsub -cwd -V -N SJrealign -S /bin/bash -q geschwind.q -l h_data=64G -pe shared 8 4_2pass.sh
#!/bin/bash

##This script now runs the 2pass star alignment using the new star index created with the splice junction file
##script is written for paired end data, with R1_001.fastq.gz and R2_001.fastq.gz for each sample

starCALL=/share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR
genomeDir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/hg19_2pass
datadir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2014-423-18765747  # this is where your fasta files are located, should contain subdirectories for each sample
outputdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar  #output directory


# Final alignments made with the index made from splice junctions "SJ.all.out.tab"
# This script goes through samples one at a time since writing sam and bam files


cd $datadir

for dir in */; do
	if [ ! -d /geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/${dir} ]; then
		mkdir ${outputdir}/${dir} 
		cd ${dir}
		for readpair in `ls *R1_001.fastq.gz`; do
		 	echo $readpair
		 	basereadpair=$(basename $readpair _R1_001.fastq.gz)
		 	echo $basereadpair
		 	cd ${outputdir}/${dir}
		 	$starCALL --runThreadN 8 --genomeDir $genomeDir --outFileNamePrefix ${outputdir}/${dir}/${basereadpair} --readFilesCommand gunzip -c --readFilesIn ${datadir}/${dir}${basereadpair}_R1_001.fastq.gz ${datadir}/${dir}${basereadpair}_R2_001.fastq.gz
		 	echo "STAR alignment complete"
		 	samtools view -bS ${basereadpair}Aligned.out.sam > ${basereadpair}Aligned.out.bam
		 	echo ".SAM now converted to .BAM"
		 	rm -f ${basereadpair}Aligned.out.sam
		 	cd $datadir
		done 	
	fi
done
