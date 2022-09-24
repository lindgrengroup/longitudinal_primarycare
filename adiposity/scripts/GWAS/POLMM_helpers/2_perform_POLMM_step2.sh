#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 08/02/2022

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 2
#$ -t 1-22 
#$ -N polmm_GWAS_step2
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

# Strata and clusters called by submission script
# chromosome from task ID
# path to bgen directory: UKB imputed files in /well/ukbb-wtchg/v3/imputation/

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/POLMM_helpers/2_POLMM_step2_wrapper.R \
--strata=${STRATA} \
--chr=${SGE_TASK_ID} \
--bgenDir=/well/ukbb-wtchg/v3/imputation/ \
--sampleFile=/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--step1ResultsDir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/ \
--step2ResultsDir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/POLMM_results/ 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
