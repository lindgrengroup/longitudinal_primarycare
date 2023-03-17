#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 24/08/22

#SBATCH -A lindgren.prj
#SBATCH -p long
#SBATCH -c 8
#SBATCH -J softclust_BOLT_GWAS_no_baseline_adj
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/BOLT_logs/softclust_BOLT_GWAS_no_baseline_adj-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

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
--statsFile=BOLT_results/${STRATA}_${CLUST}_no_baseline_adjustment_assoc_cal.stats.gz \
--bgenFile=${UKB_PATH}/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=${UKB_PATH}/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=BOLT_results/${STRATA}_${CLUST}_no_baseline_adjustment_assoc_imp.stats.gz \
--verboseStats

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
