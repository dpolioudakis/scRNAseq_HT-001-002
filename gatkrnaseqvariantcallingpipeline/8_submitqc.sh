#!/bin/bash

## To submit this script: qsub -cwd -S /bin/bash -V -N Submit_QC -q geschwind.q -l h_data=4G,h_rt=1:00:00 8_submitqc.sh
## code to submiut 9_QC.sh

CODEDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code
WORKDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/MergeRound1
GROUP=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/rerun_qc1.txt  ##can split up samples in multiple groups if lots of samples

cd $CODEDIR
mkdir -p $WORKDIR/log4

while read names;do

processdir=$WORKDIR/$names
echo $processdir

qsub -cwd -o $WORKDIR/log4 -e $WORKDIR/log4 -S /bin/bash -V -N QC_$names -q geschwind.q -l h_data=36G $CODEDIR/9_QC.sh $processdir $names


done<$GROUP
