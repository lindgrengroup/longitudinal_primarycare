#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 07/02/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 12
#$ -N saige_step1
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/SAIGE_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/SAIGE_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module purge
module load Anaconda3/2020.07
# source activate /well/lindgren/users/mmq446/conda/skylake/envs/RSAIGE

source activate /well/lindgren/users/mmq446/conda/skylake/envs/saige

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

mkdir /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/SAIGE_results/${STRATA}

# Strata and clusters called by submission script
# PLINK files: in /well/lindgren/UKBIOBANK/ferreira/IMPUTED_association_analysis/grm_files/
# Phenotype files: in /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/traits_for_GWAS/
# Separate pheno column for each cluster vs the rest (k1:k5)
# Don't add sex as covariate here because we will add that in the R wrapper
# Adjust for baseline age and age-squared

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/SAIGE_helpers/1_SAIGE_step1_wrapper.R \
--strata=${STRATA} \
--cluster=${KI} \
--adjBaseline=FALSE \
--plinkFile=/well/lindgren/UKBIOBANK/ferreira/IMPUTED_association_analysis/grm_files/ukb_cal_v2_qced_pruned \
--covars=baseline_age,age_sq,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
--sampleIDColinphenoFile=eid \
--traitType=binary \
--outputdir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/SAIGE_results/ 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
