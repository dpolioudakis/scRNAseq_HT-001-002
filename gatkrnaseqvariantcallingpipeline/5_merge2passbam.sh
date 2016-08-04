## To submit this script: qsub -cwd -S /bin/bash -V -N merge -q geschwind.q -l h_data=64G 5_merge2passbam.sh
#!/bin/bash

##This script merged 2pass bam files per sample
##Might need to update "INPUT" in line 23 for how many bam files exist per sample, bamlist will find all samples with ${sampleid}*Aligned.out.bam and then merge 

CODEDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code
WORKDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar
OUTDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam
GROUP=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/sampleid.txt   ## file with sample names
PIC=/share/apps/picard-tools-1.128/dist/picard.jar

cd $WORKDIR

while read sampleid;do

	mkdir $OUTDIR/$sampleid
	if [ ! -e $OUTDIR/${sampleid}/${sampleid}.merge.bam ]; then
	##Find all the deduped files with the sample id
	bamlist=`ls -l $(find -name "${sampleid}*Aligned.out.bam") | awk '{print $9}' | tr '\n' '\t'`
	bams=(${bamlist//'\t'/})
	##Merge
	java -Xmx2g -Djava.io.tmpdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/tmp -jar ${PIC} MergeSamFiles INPUT=${bams[0]} INPUT=${bams[1]} INPUT=${bams[2]} INPUT=${bams[3]} OUTPUT=$OUTDIR/${sampleid}/${sampleid}.merge.bam
	fi
done<$GROUP
		  