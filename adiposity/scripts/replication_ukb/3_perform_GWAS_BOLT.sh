#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q long.qc
#$ -pe shmem 8
#$ -t 1-16 -tc 4 
#$ -N GWAS_hormones
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

/apps/well/bolt-lmm/2.3.2/bolt \
--bed=/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chr{1:22}_v2.bed \
--bim=/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chr{1:22}_v2.bim \
--fam=/well/lindgren/UKBIOBANK/samvida/general_resources/BOLT-LMM/ukb11867_cal_chr1_v2_s488363_COL6NUMERIC.fam \
--remove=/well/lindgren/UKBIOBANK/samvida/general_resources/BOLT-LMM/bolt.in_plink_but_not_imputed_AUTOSOMES.FID_IID.968.txt \
--phenoFile=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/traits_for_gwas/qcd_${STRATA_NAME}_cross_sectional.txt \
--phenoCol=adj_trait \
--covarFile=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/sample_qc/${STRATA_NAME}_ids_passed_qc.txt \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_results/${STRATA_NAME}_assoc_cal.stats.gz \
--bgenFile=/well/lindgren/UKBIOBANK/DATA/IMPUTATION/ukb_imp_chr{1:22}_v3.bgen \
--sampleFile=/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--statsFileBgenSnps=/well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_results/${STRATA_NAME}_assoc_imp.stats.gz \
--verboseStats

cd /well/lindgren/UKBIOBANK/samvida/hormone_ehr/GWAS/BOLT_results
cat ./${STRATA_NAME}_assoc_cal.stats.gz \
./${STRATA_NAME}_assoc_imp.stats.gz \
> ./${STRATA_NAME}_assoc.stats.gz

# Create log file for filtering results
LOG_FILE="../log_files/${STRATA_NAME}.txt"
rm $LOG_FILE
touch $LOG_FILE

# Create temporary directory for filtering
rm -r ./${STRATA_NAME}_tmp_QC
mkdir ./${STRATA_NAME}_tmp_QC
gunzip ./${STRATA_NAME}_assoc.stats.gz 

PREQC=$(wc -l < ./${STRATA_NAME}_assoc.stats)
printf "** Phenotype and Strata: ${STRATA_NAME}\n" >> $LOG_FILE 
printf "\t # SNPs pre-QC: $((${PREQC}-1)) \n" >> $LOG_FILE

# Location of filtering files
FILTER_LOC="/well/lindgren/UKBIOBANK/samvida/full_primary_care/GWAS/snps_passed_QC_211209/"

# MAF filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_mfi_maf_QC.txt ./${STRATA_NAME}_assoc.stats \
> ./${STRATA_NAME}_tmp_QC/tmp_passed_maf.txt

TMPQC=$(wc -l < ./${STRATA_NAME}_tmp_QC/tmp_passed_maf.txt)
printf "\t # SNPs passed MAF > 0.01: $((${TMPQC}-1)) \n" >> $LOG_FILE

# INFO filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_info_QC.txt \
./${STRATA_NAME}_tmp_QC/tmp_passed_maf.txt \
> ./${STRATA_NAME}_tmp_QC/tmp_passed_info.txt

TMPQC=$(wc -l < ./${STRATA_NAME}_tmp_QC/tmp_passed_info.txt)
printf "\t # SNPs passed INFO > 0.8: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Genotyping filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_geno_QC.txt \
./${STRATA_NAME}_tmp_QC/tmp_passed_info.txt \
> ./${STRATA_NAME}_tmp_QC/tmp_passed_geno.txt

TMPQC=$(wc -l < ./${STRATA_NAME}_tmp_QC/tmp_passed_geno.txt)
printf "\t # SNPs passed missingness < 0.05: $((${TMPQC}-1)) \n" >> $LOG_FILE

# HWE filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_hwe_QC.txt \
./${STRATA_NAME}_tmp_QC/tmp_passed_geno.txt \
> ./${STRATA_NAME}_tmp_QC/tmp_passed_hwe.txt

TMPQC=$(wc -l < ./${STRATA_NAME}_tmp_QC/tmp_passed_hwe.txt)
printf "\t # SNPs passed HWE pval > 1E-06: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Biallelic filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_biallelic_QC.txt \
./${STRATA_NAME}_tmp_QC/tmp_passed_hwe.txt \
> ./${STRATA_NAME}_tmp_QC/tmp_passed_biallelic.txt

TMPQC=$(wc -l < ./${STRATA_NAME}_tmp_QC/tmp_passed_biallelic.txt)
printf "\t # SNPs passed bi-allelic QC: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Save results
cp ./${STRATA_NAME}_tmp_QC/tmp_passed_biallelic.txt ./${STRATA_NAME}_filtered.txt
gzip ./${STRATA_NAME}_filtered.txt
rm -r ./${STRATA_NAME}_tmp_QC

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
