#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 26/08/22

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -N h2_rg_calc
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

/apps/well/bolt-lmm/2.3.2/bolt \
--bed=/well/lindgren-ukbb/projects/ukbb-11867/DATA/CALLS/ukb_cal_chr{1:22}_v2.bed \
--bim=/well/lindgren-ukbb/projects/ukbb-11867/DATA/CALLS/ukb_snp_chr{1:22}_v2.bim \
--fam=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/BOLT-LMM/ukb11867_cal_chr1_v2_s488363_COL6NUMERIC.fam \
--remove=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/BOLT-LMM/bolt.in_plink_but_not_imputed_AUTOSOMES.FID_IID.968.txt \
--phenoFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_traits_${STRATA}.txt \
--phenoCol=lmm_intercepts \
--phenoCol=lmm_slopes_adj_int \
--covarFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_GWAS/lmm_traits_${STRATA}.txt \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--reml \
--numThreads=8 \
--verboseStats \
2>&1 | tee /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/${STRATA}/h2rg_BOLTREML.log

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
