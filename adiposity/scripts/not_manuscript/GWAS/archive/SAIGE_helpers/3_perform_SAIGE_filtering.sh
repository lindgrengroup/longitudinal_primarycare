#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 08/02/2022

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N filter_GWAS
#$ -t 2-5 -tc 1 
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_logs/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_logs/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA="Weight_sex_comb"
KI=${SGE_TASK_ID}

# Location of filtering files
FILTER_LOC="/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/snps_passed_QC_211209/"

# Create log file for filtering results
LOG_FILE="/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_logs/filtering_${STRATA}_k${KI}.txt"
rm $LOG_FILE
touch $LOG_FILE

# Create temporary directory for filtering
cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/${STRATA}
rm -r k${KI}_tmp_QC
mkdir k${KI}_tmp_QC

# Go through each chromosome and convert chr:pos_A1_A2 varids to rsids
for CHR in {1..22}
do
	sort -k3 k${KI}_chr${CHR}_gwas_results.txt > k${KI}_tmp_QC/tmp_gwas_sorted.txt
	sort -k1 /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_mfi_chr${CHR}_v3.txt > k${KI}_tmp_QC/tmp_mfi_sorted.txt

	join -1 3 -2 1 -o 2.2,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21 \
	k${KI}_tmp_QC/tmp_gwas_sorted.txt k${KI}_tmp_QC/tmp_mfi_sorted.txt \
	> k${KI}_tmp_QC/tmp_chr${CHR}_gwas_results.txt
done

# Merge results from all chromosomes for relevant strata and cluster
# Write header
head -1 k${KI}_chr1_gwas_results.txt > k${KI}_tmp_QC/tmp_gwas_results.txt
sed -i '1s/^/rsid /' k${KI}_tmp_QC/tmp_gwas_results.txt
cat k${KI}_tmp_QC/tmp_chr*_gwas_results.txt >> k${KI}_tmp_QC/tmp_gwas_results.txt

PREQC=$(wc -l < k${KI}_tmp_QC/tmp_gwas_results.txt)
printf "** Strata and Cluster: ${STRATA}, k${KI}\n" >> $LOG_FILE 
printf "\t # SNPs pre-QC: $((${PREQC}-1)) \n" >> $LOG_FILE

# MAF filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_mfi_maf_QC.txt \
k${KI}_tmp_QC/tmp_gwas_results.txt \
> k${KI}_tmp_QC/tmp_passed_maf.txt

TMPQC=$(wc -l < k${KI}_tmp_QC/tmp_passed_maf.txt)
printf "\t # SNPs passed MAF > 0.01: $((${TMPQC}-1)) \n" >> $LOG_FILE

# INFO filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_info_QC.txt \
k${KI}_tmp_QC/tmp_passed_maf.txt \
> k${KI}_tmp_QC/tmp_passed_info.txt

TMPQC=$(wc -l < k${KI}_tmp_QC/tmp_passed_info.txt)
printf "\t # SNPs passed INFO > 0.8: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Genotyping filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_geno_QC.txt \
k${KI}_tmp_QC/tmp_passed_info.txt \
> k${KI}_tmp_QC/tmp_passed_geno.txt

TMPQC=$(wc -l < k${KI}_tmp_QC/tmp_passed_geno.txt)
printf "\t # SNPs passed missingness < 0.05: $((${TMPQC}-1)) \n" >> $LOG_FILE

# HWE filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_hwe_QC.txt \
k${KI}_tmp_QC/tmp_passed_geno.txt \
> k${KI}_tmp_QC/tmp_passed_hwe.txt

TMPQC=$(wc -l < k${KI}_tmp_QC/tmp_passed_hwe.txt)
printf "\t # SNPs passed HWE pval > 1E-06: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Biallelic filter
awk -F ' ' 'NR==FNR{a[$1]; next} FNR==1 || $1 in a' \
${FILTER_LOC}/passed_biallelic_QC.txt \
k${KI}_tmp_QC/tmp_passed_hwe.txt \
> k${KI}_tmp_QC/tmp_passed_biallelic.txt

TMPQC=$(wc -l < k${KI}_tmp_QC/tmp_passed_biallelic.txt)
printf "\t # SNPs passed bi-allelic QC: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Save results
rm ${STRATA}_k${KI}_gwas_results.txt.gz
cp k${KI}_tmp_QC/tmp_passed_biallelic.txt ${STRATA}_k${KI}_gwas_results.txt
gzip ${STRATA}_k${KI}_gwas_results.txt
rm -r k${KI}_tmp_QC

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

# Strata and clusters called by submission script
# path to where plots should be stored
mkdir /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/plots/${STRATA}/k${KI}

Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/SAIGE_helpers/3_SAIGE_filtering_wrapper.R \
--inputFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/${STRATA}/${STRATA}_k${KI}_gwas_results.txt.gz \
--logFile=$LOG_FILE \
--outputFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/SAIGE_results/${STRATA}/${STRATA}_k${KI}_final.txt \
--outPlotDir=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/plots/${STRATA}/k${KI}/ 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
