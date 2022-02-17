#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N filter_BOLT_results
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA}
rm ${STRATA}_${PARAMETER}_assoc.stats.gz
cat ${STRATA}_${PARAMETER}_assoc_cal.stats.gz \
${STRATA}_${PARAMETER}_assoc_imp.stats.gz \
> ${STRATA}_${PARAMETER}_assoc.stats.gz

# Create log file for filtering results
LOG_FILE="/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/log_files/post_GWAS/${STRATA}_${PARAMETER}.txt"
rm $LOG_FILE
touch $LOG_FILE

# Create temporary directory for filtering
rm -r ./${STRATA}_${PARAMETER}_tmp_QC
mkdir ./${STRATA}_${PARAMETER}_tmp_QC
gunzip ${STRATA}_${PARAMETER}_assoc.stats.gz 

PREQC=$(wc -l < ${STRATA}_${PARAMETER}_assoc.stats)
printf "** Phenotype and Strata: ${STRATA}_${PARAMETER}\n" >> $LOG_FILE 
printf "\t # SNPs pre-QC: $((${PREQC}-1)) \n" >> $LOG_FILE

# Location of filtering files
FILTER_LOC="/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/snps_passed_QC_211209/"

# MAF filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_mfi_maf_QC.txt ${STRATA}_${PARAMETER}_assoc.stats \
> ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_maf.txt

TMPQC=$(wc -l < ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_maf.txt)
printf "\t # SNPs passed MAF > 0.01: $((${TMPQC}-1)) \n" >> $LOG_FILE

# INFO filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_info_QC.txt \
./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_maf.txt \
> ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_info.txt

TMPQC=$(wc -l < ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_info.txt)
printf "\t # SNPs passed INFO > 0.8: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Genotyping filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_geno_QC.txt \
./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_info.txt \
> ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_geno.txt

TMPQC=$(wc -l < ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_geno.txt)
printf "\t # SNPs passed missingness < 0.05: $((${TMPQC}-1)) \n" >> $LOG_FILE

# HWE filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_hwe_QC.txt \
./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_geno.txt \
> ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_hwe.txt

TMPQC=$(wc -l < ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_hwe.txt)
printf "\t # SNPs passed HWE pval > 1E-06: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Biallelic filter
awk -F '\t' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_biallelic_QC.txt \
./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_hwe.txt \
> ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_biallelic.txt

TMPQC=$(wc -l < ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_biallelic.txt)
printf "\t # SNPs passed bi-allelic QC: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Save results
rm ${STRATA}_${PARAMETER}_filtered.txt.gz
cp ./${STRATA}_${PARAMETER}_tmp_QC/tmp_passed_biallelic.txt ${STRATA}_${PARAMETER}_filtered.txt
gzip ${STRATA}_${PARAMETER}_filtered.txt
rm -r ./${STRATA}_${PARAMETER}_tmp_QC

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

rm -r /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/plots/${STRATA}/${PARAMETER}/ 
mkdir /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/plots/${STRATA}/${PARAMETER}/ 

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_helpers/2_BOLT_filtering_wrapper.R \
--inputFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA}/${STRATA}_${PARAMETER}_filtered.txt.gz \
--logFile=$LOG_FILE \
--outputFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA}/${STRATA}_${PARAMETER}_final.txt \
--outPlotDir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/plots/${STRATA}/${PARAMETER}/ 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
