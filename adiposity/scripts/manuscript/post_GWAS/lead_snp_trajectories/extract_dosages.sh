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

UKB_PATH="" # REDACTED 

/apps/well/qctool/2.0.1/qctool \
-g ${UKB_PATH}/IMPUTATION/ukb_imp_chr${CHR}_v3.bgen \
-s ${UKB_PATH}/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
-condition-on ${VARID} \
-os ${VARID}_dosages.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

