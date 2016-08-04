## Script to submit an Rscript to Orion queue
## 2_SJ.out.ALL.R is used to created a combined splice junction file for 2pass star alignment
## Use this qsub: qsub -cwd -V -N RunR -S /bin/bash -l h_data=16G 1_RunSJ.sh

#!/bin/bash



cd /geschwindlabshares/eQTL/data/FetalBrainRNASeq/2passStar/code/   #code directory where 2_SJ.out.ALL.R is located

Rscript --vanilla 2_SJ.out.ALL.R