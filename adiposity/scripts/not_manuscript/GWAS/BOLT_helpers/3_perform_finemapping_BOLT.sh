#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/12/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N finemap
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/BOLT_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/BOLT_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

module purge
cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline
source activate pipeline

pipeline.r --outname=${STRATA}_${PARAMETER} \
--outdir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS/${STRATA} \
--input-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/0_input_jsons/${STRATA}_inputs.json \
--job-dependencies-json=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/pipeline/json/job_dependencies.json \
--sumstats=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/BOLT_results/${STRATA}_${PARAMETER}_final.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
