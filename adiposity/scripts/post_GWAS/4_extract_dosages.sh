#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 28/02/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N extract_dosages

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

/apps/well/qctool/2.0.1/qctool \
-g /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${CHR}_v3.bgen \
-s /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
-condition-on ${VARID} \
-os /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sample_variant_counts/${VARID}_dosages.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

