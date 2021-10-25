#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -t 1-24 -tc 8 
#$ -N spline_coeffs_BOLT_GWAS
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`
FILE_NAME="$(echo $STRATA_NAME | cut -d'.' -f1)"
TERM_NAME="$(echo $STRATA_NAME | cut -d'.' -f2)"

/apps/well/bolt-lmm/2.3.2/bolt \
--bed=/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr{1:22}_v2.bed \
--bim=/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr{1:22}_v2.bim \
--fam=/well/lindgren/UKBIOBANK/samvida/general_resources/BOLT-LMM/ukb11867_cal_chr1_v2_s488363_COL6NUMERIC.fam \
--remove=/well/lindgren/UKBIOBANK/samvida/general_resources/BOLT-LMM/bolt.in_plink_but_not_imputed_AUTOSOMES.FID_IID.968.txt \
--phenoFile=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/traits_for_gwas/${FILE_NAME}.txt \
--phenoCol=${TERM_NAME} \
--covarFile=ids_passed_qc_211015.txt \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA_NAME}_assoc_cal.stats.gz \
--bgenFile=/well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF=0.001 \
--bgenMinINFO=0.8 \
--maxMissingPerSnp=0.05 \
--maxMissingPerIndiv=0.05 \
--sampleFile=/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA_NAME}_assoc_imp.stats.gz \
--verboseStats

cat BOLT_results/${STRATA_NAME}_assoc_cal.stats.gz BOLT_results/${STRATA_NAME}_assoc_imp.stats.gz > BOLT_results/${STRATA_NAME}_assoc.stats.gz