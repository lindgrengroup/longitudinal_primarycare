#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -N lmm_BOLT_GWAS
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_logs/
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
--phenoFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/traits_for_gwas/${PARAMETER}_${STRATA}.txt \
--phenoCol=adj_trait \
--covarFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA}_ids_passed_qc.txt \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/${STRATA}_${PARAMETER}_assoc_cal.stats.gz \
--bgenFile=/well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/BOLT_results/${STRATA}_${PARAMETER}_assoc_imp.stats.gz \
--verboseStats

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
