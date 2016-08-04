#!/bin/bash

## Script to submit an Rscript to Orion queue
## 2_SJ.out.ALL.R is used to created a combined splice junction file for 2pass star alignment
## Use this qsub:
# qsub -cwd -V -N varGATK1 -S /bin/bash -l h_data=16G -o logs/Variant_GATK_1_RunSJ_QSUB_$(date +%Y%m%d).log -e logs/Variant_GATK_1_RunSJ_QSUB_$(date +%Y%m%d).error Variant_GATK_1_RunSJ.sh

/share/apps/R-3.2.2/bin/Rscript --vanilla Variant_GATK_2_SJ.out.ALL.R
