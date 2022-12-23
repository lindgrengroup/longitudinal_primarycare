#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N finemap
#$ -o BOLT_logs/
#$ -e BOLT_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

GEN_RESOURCES_PATH="" # REDACTED

module purge
source activate pipeline

pipeline.r --outname=${STRATA}_${PARAMETER} \
--outdir=post_GWAS/${STRATA} \
--input-json=post_GWAS/0_input_jsons/${STRATA}_inputs.json \
--job-dependencies-json=${GEN_RESOURCES_PATH}/pipeline/json/job_dependencies.json \
--sumstats=BOLT_results/${STRATA}_${PARAMETER}_final.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
