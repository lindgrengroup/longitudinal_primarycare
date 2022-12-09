#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -t 1-6 
#$ -N GWAS_adipo_cross_sec_no_gp
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

/apps/well/bolt-lmm/2.3.2/bolt \
--bed=/well/lindgren-ukbb/projects/ukbb-11867/DATA/CALLS/ukb_cal_chr{1:22}_v2.bed \
--bim=/well/lindgren-ukbb/projects/ukbb-11867/DATA/CALLS/ukb_snp_chr{1:22}_v2.bim \
--fam=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/BOLT-LMM/ukb11867_cal_chr1_v2_s488363_COL6NUMERIC.fam \
--remove=/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/BOLT-LMM/bolt.in_plink_but_not_imputed_AUTOSOMES.FID_IID.968.txt \
--phenoFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/traits_for_gwas/qcd_${STRATA_NAME}_cross_sectional.txt \
--phenoCol=adj_trait \
--covarFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/sample_qc/${STRATA_NAME}_ids_passed_qc.txt \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/BOLT_results/${STRATA_NAME}_assoc_cal.stats.gz \
--bgenFile=/well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/BOLT_results/${STRATA_NAME}_assoc_imp.stats.gz \
--verboseStats

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/GWAS/BOLT_results
cat ./${STRATA_NAME}_assoc_cal.stats.gz \
./${STRATA_NAME}_assoc_imp.stats.gz \
> ./${STRATA_NAME}_assoc.stats.gz

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
