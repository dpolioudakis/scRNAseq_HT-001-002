## use this to submit 7_2pass_2.sh
## To submit this script: qsub -cwd -S /bin/bash -V -N Submit_var -q geschwind.q -l h_data=4G,h_rt=3:00:00 6_sub2pass_2.sh
#!/bin/bash

##This script lets samples run in groups and will hold samples until group before is finished running
##Need to create group.txt which lists files of groups ie group1.txt,group2.txt etc
##THen you need to have each group1.txt list 10-15 sample names


CODEDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code
WORKDIR=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/Merged2passBam
GROUP=/geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/group.txt
JOBNAME=var2pass 				## a unique identifier for this job submission


cd $CODEDIR

mkdir -p $WORKDIR/log

while read groups;do

	num=$( echo "$groups" | cut -c 68-69 )
	echo $num

	cd $WORKDIR

	if [ "$groups" = "batch01.txt" ]; then

		while read names;do

		processdir=$WORKDIR/$names

		mkdir -p $processdir

		echo $processdir
		
		qsub -cwd -o $WORKDIR/log -e $WORKDIR/log -V -S /bin/bash -q geschwind.q -N var_${names}_${JOBNAME}_${num} -l h_data=36G,h_rt=24:00:00 ${CODEDIR}/7_2pass_2.sh $processdir $names

		done<$groups

	else

		while read names;do

                processdir=$WORKDIR/$names

                mkdir -p $processdir

                echo $processdir

		hold_id=$((10#$num -1))
		
                qsub -hold_jid *${JOBNAME}*${hold_id}  -cwd -o $WORKDIR/log -e $WORKDIR/log -V -S /bin/bash -q geschwind.q -N var_${names}_${JOBNAME}_${num} -l h_data=36G,h_rt=24:00:00 ${CODEDIR}/7_2pass_2.sh $processdir $names

                done<$groups

	fi


done<$GROUP
