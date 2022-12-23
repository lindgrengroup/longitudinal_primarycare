#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 24/08/22

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -N soft_clust_BOLT_GWAS
#$ -o BOLT_logs
#$ -e BOLT_logs
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

UKB_PATH="" # REDACTED
GEN_RESOURCES_PATH="" # REDACTED

mkdir BOLT_results/${STRATA}

/apps/well/bolt-lmm/2.3.2/bolt \
--bed=${UKB_PATH}/CALLS/ukb_cal_chr{1:22}_v2.bed \
--bim=${UKB_PATH}/CALLS/ukb_snp_chr{1:22}_v2.bim \
--fam=${GEN_RESOURCES_PATH}/BOLT-LMM/ukb11867_cal_chr1_v2_s488363_COL6NUMERIC.fam \
--remove=${GEN_RESOURCES_PATH}/BOLT-LMM/bolt.in_plink_but_not_imputed_AUTOSOMES.FID_IID.968.txt \
--phenoFile=traits_for_GWAS/${STRATA}_soft_clust_probs.txt \
--phenoCol=${CLUST} \
--covarFile=sample_qc/${STRATA}_covariates.txt \
--qCovarCol=baseline_trait \
--qCovarCol=baseline_age \
--qCovarCol=age_sq \
--qCovarCol=FUyrs \
--qCovarCol=FU_n \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--covarCol=UKB_assmt_centre \
--covarCol=sex \
--covarMaxLevels=30 \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=BOLT_results/${STRATA}_${CLUST}_assoc_cal.stats.gz \
--bgenFile=${UKB_PATH}/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=${UKB_PATH}/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=BOLT_results/${STRATA}_${CLUST}_assoc_imp.stats.gz \
--verboseStats

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
