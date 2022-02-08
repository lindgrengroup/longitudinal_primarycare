#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 08/02/2022

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 2
#$ -t 1-22 
#$ -N cluster_GWAS_step2
#$ -o /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_logs/
#$ -e /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge
module load Anaconda3/2020.07
source activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

# Strata and clusters called by submission script
# chromosome from task ID
# path to bgen directory: UKB imputed files in /well/ukbb-wtchg/v3/imputation/
# path to sample QC ids: /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/
# path to results from step 1: /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/
# path to where output should be stored: /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/

Rscript /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/SAIGE_step2_wrapper.R \
--strata=${STRATA} \
--cluster=k${KI} \
--chr=${SGE_TASK_ID} \
--bgenDir=/well/ukbb-wtchg/v3/imputation/ \
--sampleDir=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/ \
--step1ResultsDir=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/ \
--step2ResultsDir=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/SAIGE_results/ 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
