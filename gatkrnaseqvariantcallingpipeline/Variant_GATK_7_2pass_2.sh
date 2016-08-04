#!/bin/bash

##Define Input Variables and Functions
PROCESSDIR=$1 	## /geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/{dir} 
NAME=$2	## directory name

##This script adds read groups, marks duplicates, sorts and reorders the bam file
##Then runs SplitNCigarReads which will split reads wiht N's in the cigar strings and creates k+1 new reads (where k is the number of N cigar elements) that correspond to the segments of the original read beside/between the splicing events represented by the Ns in the original CIGAR
##Then base recalibration is performed
##Finally variant calling and variant filtering

PIC=/share/apps/picard-tools-1.128/dist/picard.jar
gatk=/home/rwalker/bin/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar
outputdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam
hg19=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/Genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.ercc.fa
snp=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/Genome/dbsnp_GRCh37p19_b146common_all.vcf




cd $PROCESSDIR

## Add read groups, sort, mark duplicates, and create index Using Picard
java -jar -Xmx4g -Djava.io.tmpdir=${outputdir}/tmp $PIC AddOrReplaceReadGroups I=${NAME}.merge.bam O=${NAME}_rg_added_sorted.bam SO=coordinate RGID=1 RGLB=library RGPL=illumina RGPU=machine RGSM=sample 
java -jar -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp $PIC MarkDuplicates I=${NAME}_rg_added_sorted.bam O=${NAME}_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${NAME}_output.metrics 
## Sort bam file
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $PIC SortSam INPUT=${NAME}_dedupped.bam OUTPUT=${NAME}_sorted.bam SORT_ORDER=coordinate TMP_DIR=${outputdir}/tmp
## Reorder bam file
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $PIC ReorderSam INPUT=${NAME}_sorted.bam OUTPUT=${NAME}_reordered.bam REFERENCE=$hg19 TMP_DIR=${outputdir}/tmp
## index bam file
samtools index ${NAME}_reordered.bam
## GATK Split'N'Trim and reassign mapping qualities
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $gatk -T SplitNCigarReads -R $hg19 -I ${NAME}_reordered.bam -o ${NAME}_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 
## Base Recalibration
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $gatk -T BaseRecalibrator -R $hg19 -I ${NAME}_split.bam  -knownSites $snp -o ${NAME}_recal_data.table 
## Recalibrate Bam
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $gatk -T PrintReads -R $hg19 -I ${NAME}_split.bam -BQSR ${NAME}_recal_data.table -o ${NAME}_recalibration.bam
## Variant Calling
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $gatk -T HaplotypeCaller -R $hg19 -I ${NAME}_recalibration.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o ${NAME}.vcf 
## Variant Filtering
java -Xmx16g -Djava.io.tmpdir=${outputdir}/tmp -jar $gatk -T VariantFiltration -R $hg19 -V ${NAME}.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o ${NAME}_filtered.vcf 
##remove intermediate bam files
rm ${NAME}_rg_added_sorted.bam
rm ${NAME}_dedupped.bam
rm ${NAME}_dedupped.bai
rm ${NAME}_sorted.bam
rm ${NAME}_reordered.bam
rm ${NAME}_reordered.bam.bai
rm ${NAME}_split.bam
rm ${NAME}_split.bai
rm ${NAME}_recal_data.table
rm ${NAME}.vcf
rm ${NAME}.vcf.idx






