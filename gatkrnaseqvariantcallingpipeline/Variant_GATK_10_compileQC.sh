#!/bin/bash

## To submit: qsub -cwd -V -N picard -S /bin/bash -q geschwind.q -l h_data=16G,h_rt=3:00:00 10_compileQC.sh


## Set location for output and variables
codedir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code
workingdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam
sampid=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/sampleid.txt
qcdir=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam/QC


mdkir $qcdir

echo "FileName CATEGORY TOTAL_READS	PF_READS PCT_PF_READS PF_NOISE_READS PF_READS_ALIGNED PCT_PF_READS_ALIGNED PF_ALIGNED_BASES PF_HQ_ALIGNED_READS PF_HQ_ALIGNED_BASES PF_HQ_ALIGNED_Q20_BASES PF_HQ_MEDIAN_MISMATCHES PF_MISMATCH_RATE PF_HQ_ERROR_RATE PF_INDEL_RATE MEAN_READ_LENGTH READS_ALIGNED_IN_PAIRS PCT_READS_ALIGNED_IN_PAIRS BAD_CYCLES STRAND_BALANCE PCT_CHIMERAS PCT_ADAPTER SAMPLE LIBRARY READ_GROUP" > ${qcdir}/AlignmentSummary.txt
echo "FileName READ_END CYCLE PCT_A PCT_C PCT_G PCT_T PCT_N" > ${qcdir}/Base_distribution.txt
echo "GC WINDOWS READ_STARTS MEAN_BASE_QUALITY NORMALIZED_COVERAGE ERROR_BAR_WIDTH" > ${qcdir}/RNAseqGC.txt
echo "WINDOW_SIZE TOTAL_CLUSTERS ALIGNED_READS AT_DROPOUT GC_DROPOUT" > ${qcdir}/GCsummary.txt
echo "FileName MEDIAN_INSERT_SIZE MEDIAN_ABSOLUTE_DEVIATION MIN_INSERT_SIZE MAX_INSERT_SIZE MEAN_INSERT_SIZE STANDARD_DEVIATION READ_PAIRS PAIR_ORIENTATION WIDTH_OF_10_PERCENT WIDTH_OF_20_PERCENT WIDTH_OF_30_PERCENT WIDTH_OF_40_PERCENT WIDTH_OF_50_PERCENT WIDTH_OF_60_PERCENT WIDTH_OF_70_PERCENT WIDTH_OF_80_PERCENT WIDTH_OF_90_PERCENT WIDTH_OF_99_PERCENT SAMPLE LIBRARY READ_GROUP" >${qcdir}/InsertSizeSummary.txt 
echo "FileName PF_BASES PF_ALIGNED_BASES RIBOSOMAL_BASES CODING_BASES UTR_BASES INTRONIC_BASES INTERGENIC_BASES IGNORED_READS CORRECT_STRAND_READS INCORRECT_STRAND_READS PCT_RIBOSOMAL_BASES PCT_CODING_BASES PCT_UTR_BASES PCT_INTRONIC_BASES PCT_INTERGENIC_BASES PCT_MRNA_BASES PCT_USABLE_BASES PCT_CORRECT_STRAND_READS MEDIAN_CV_COVERAGE MEDIAN_5PRIME_BIAS MEDIAN_3PRIME_BIAS MEDIAN_5PRIME_TO_3PRIME_BIAS SAMPLE LIBRARY READ_GROUP" > ${qcdir}/RNAseqQC.txt
echo "" > ${qcdir}/TranscriptCoverage.txt
echo "FileName UniquelyMappedReads" > ${qcdir}/UniquelyMappedReads.txt


cd $workingdir

while read samp;do
	cd $workingdir/${samp}
	
 	echo "Getting alignment summary"
 	var=`sed -n 10p ${workingdir}/${samp}/${samp}_alignment_summary.txt`
 	echo ${samp} ${var} >> ${qcdir}/AlignmentSummary.txt

 	echo "Getting base distribution summary"
 	var=`sed -n 8,107p ${workingdir}/${samp}/${samp}_base_distribution.txt`
 	echo ${samp} ${var} >> ${qcdir}/Base_distribution.txt 	
	
	echo "Getting gcbias stats"
	var=`sed -n 8,1088p ${workingdir}/${samp}/${samp}_gcbias_stats.txt`
	echo ${samp} ${var} >> ${qcdir}/RNAseqGC.txt

	echo "Getting gcbias summary"
	var=`sed -n 8p ${workingdir}/${samp}/${samp}_gcbias_summary.txt`
	echo ${samp} ${var} >> ${qcdir}/GCsummary.txt

    echo "Getting insert size stats"
    var=`sed -n 8p ${workingdir}/${samp}/${samp}_insert_size.txt`
    echo ${samp} ${var} >> ${qcdir}/InsertSizeSummary.txt

	echo "Getting rnaseq stats"
	var=`sed -n 8p ${workingdir}/${samp}/${samp}_rnaseq_stats.txt`
	echo ${samp} ${var} >> ${qcdir}/RNAseqQC.txt
	
	var=`sed -n 11,112p ${workingdir}/${samp}/${samp}_rnaseq_stats.txt`
	echo ${samp} ${var} >> ${qcdir}/TranscriptCoverage.txt
	
	echo "Getting uniquely mapped reads"
	var=`sed -n 9p ${workingdir}/${samp}/${samp}Log.final.out`
	echo ${samp} ${var} >> ${qcdir}/UniquelyMappedReads.txt
	
	cd $workingdir

done<$sampid















