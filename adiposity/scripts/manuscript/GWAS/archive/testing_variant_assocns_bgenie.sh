#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/11/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -N metabo_snp_assocn_bgenie
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load BGEN/1.1.6-GCCcore-8.3.0

bgenie \
--bgen=/well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr21_v3.bgen \
--pheno=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/all_traits_for_bgenie.txt \
--covar=/well/lindgren/UKBIOBANK/samvida/general_resources/genetic_covariates_for_bgenie.txt \
--out=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/BGENIE_results/chr21 \
--pvals \
--include_rsids=/well/lindgren/UKBIOBANK/samvida/full_primary_care/GWAS/metabo_endo_snps_passed_QC_211102/passed_QC.txt \
--thread=8 
