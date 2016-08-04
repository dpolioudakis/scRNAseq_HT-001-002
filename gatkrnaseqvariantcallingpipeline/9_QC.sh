#!/bin/bash

##QC stats
##Define Input Variables and Functions
PROCESSDIR=$1 	## /geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/{dir} 
NAME=$2	## directory name

FASTA=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/Genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.reordered.ercc.fa
PIC=/share/apps/picard-tools-1.128/dist/picard.jar
REFFLAT=/home/rwalker/fetal_brain_eqtl/genePred.v19
outputdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam


##QC on 2passStar using Picard Tools 
##Uses default values for all programs unless stated below
##
##CollectAlignmentSummaryMetrics options: ADAPTER_SEQUENCE=null
##CollectInsertSizeMetrics options: HISTOGRAM_FILE=${PROCESSDIR}/${NAME}.insertsizeHist.pdf
##QualityScoreDistribution options: CHART_OUTPUT=${PROCESSDIR}/${NAME}.qualscore.pdf
##MeanQualityByCycle options: CHART_OUTPUT=${OUTDIR}/mean_quality_chart.pdf
##CollectBaseDistributionByCycle  options: CHART_OUTPUT=${OUTDIR}/base_distribution_chart.pdf
##CollectGcBiasMetrics options: CHART_OUTPUT=${OUTDIR}/gcbias_chart.pdf SUMMARY_OUTPUT=${OUTDIR}/gcbias_summary.txt
##CollectRnaSeqMetrics options: STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND REF_FLAT=${REFFLAT} 
##BamIndexStats options: none
##EstimateLibraryComplexity options: none


##${base}_recalibration.bam has already added read groups, marked duplicates, sorted and reordered, splitncigar, base recalibration

cd $PROCESSDIR

java -Xmx2g -jar ${PIC} CollectAlignmentSummaryMetrics 	INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_alignment_summary.txt ASSUME_SORTED=true ADAPTER_SEQUENCE=null
java -Xmx2g -jar -Djava.io.tmpdir=${outputdir}/tmp ${PIC} CollectInsertSizeMetrics 		INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_insert_size.txt ASSUME_SORTED=true HISTOGRAM_FILE=${NAME}_insertsizeHist.pdf TMP_DIR=${outputdir}/tmp
java -Xmx2g -jar -Djava.io.tmpdir=${outputdir}/tmp ${PIC} QualityScoreDistribution 		INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_quality_score.txt ASSUME_SORTED=true CHART_OUTPUT=${NAME}_quality_score_chart.pdf TMP_DIR=${outputdir}/tmp
java -Xmx2g -jar -Djava.io.tmpdir=${outputdir}/tmp ${PIC} MeanQualityByCycle 			INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_mean_quality.txt ASSUME_SORTED=true CHART_OUTPUT=${NAME}_mean_quality_chart.pdf TMP_DIR=${outputdir}/tmp
java -Xmx2g -jar -Djava.io.tmpdir=${outputdir}/tmp ${PIC} CollectBaseDistributionByCycle	INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_base_distribution.txt ASSUME_SORTED=true CHART_OUTPUT=${NAME}_base_distribution_chart.pdf TMP_DIR=${outputdir}/tmp
java -Xmx2g -jar -Djava.io.tmpdir=${outputdir}/tmp ${PIC} CollectGcBiasMetrics 			INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_gcbias_stats.txt ASSUME_SORTED=true CHART_OUTPUT=${NAME}_gcbias_chart.pdf SUMMARY_OUTPUT=${NAME}_gcbias_summary.txt TMP_DIR=${outputdir}/tmp
java -Xmx2g -jar ${PIC} CollectRnaSeqMetrics 		INPUT=${NAME}_recalibration.bam REFERENCE_SEQUENCE=${FASTA} OUTPUT=${NAME}_rnaseq_stats.txt ASSUME_SORTED=true STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND REF_FLAT=${REFFLAT}
java -Xmx2g -jar ${PIC} BamIndexStats 				INPUT=${NAME}_recalibration.bam
java -Xmx2g -jar ${PIC} EstimateLibraryComplexity	INPUT=${NAME}_recalibration.bam OUTPUT=${NAME}_lib_complexity.txt

