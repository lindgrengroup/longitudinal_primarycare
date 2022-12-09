#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 01/06/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 3
#$ -t 1-22 -tc 5
#$ -N prune_extract
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/obesity_genotypes

mkdir tmp

# Get obesity-associated SNPs on this chromosome
awk -v pat="chr${SGE_TASK_ID}" '$0 ~ pat{print $4}' \
/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed \
> tmp/chr${SGE_TASK_ID}_all_variants.txt
sort -o tmp/chr${SGE_TASK_ID}_all_variants.txt tmp/chr${SGE_TASK_ID}_all_variants.txt

# Ensure that the selected variants pass all QC filters

# Create log file for filtering results
LOG_FILE="logs/chr${SGE_TASK_ID}_filtering.txt"
rm $LOG_FILE
touch $LOG_FILE

# Create temporary directory for filtering
rm -r chr${SGE_TASK_ID}_tmp_QC
mkdir chr${SGE_TASK_ID}_tmp_QC

PREQC=$(wc -l < tmp/chr${SGE_TASK_ID}_all_variants.txt)
printf "\t # SNPs pre-QC: $((${PREQC}-1)) \n" >> $LOG_FILE

# Location of filtering files
FILTER_LOC="/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/GWAS/snps_passed_QC_220623/"

# MAF 0.001 filter
comm -12 ${FILTER_LOC}/passed_mfi_maf_QC.txt tmp/chr${SGE_TASK_ID}_all_variants.txt \
> ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_maf.txt
TMPQC=$(wc -l < ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_maf.txt)
printf "\t # SNPs passed MAF > 0.001: $((${TMPQC}-1)) \n" >> $LOG_FILE

# INFO filter
comm -12 ${FILTER_LOC}/passed_info_QC.txt chr${SGE_TASK_ID}_tmp_QC/tmp_passed_maf.txt \
> ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_info.txt
TMPQC=$(wc -l < ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_info.txt)
printf "\t # SNPs passed INFO > 0.8: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Genotyping filter
comm -12 ${FILTER_LOC}/passed_geno_QC.txt chr${SGE_TASK_ID}_tmp_QC/tmp_passed_info.txt \
> ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_geno.txt
TMPQC=$(wc -l < ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_geno.txt)
printf "\t # SNPs passed missingness < 0.05: $((${TMPQC}-1)) \n" >> $LOG_FILE

# HWE filter
comm -12 ${FILTER_LOC}/passed_hwe_QC.txt chr${SGE_TASK_ID}_tmp_QC/tmp_passed_geno.txt \
> ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_hwe.txt
TMPQC=$(wc -l < ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_hwe.txt)
printf "\t # SNPs passed HWE pval > 1E-06: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Biallelic filter
comm -12 ${FILTER_LOC}/passed_biallelic_QC.txt chr${SGE_TASK_ID}_tmp_QC/tmp_passed_hwe.txt \
> ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_biallelic.txt
TMPQC=$(wc -l < ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_biallelic.txt)
printf "\t # SNPs passed bi-allelic QC: $((${TMPQC}-1)) \n" >> $LOG_FILE

# Save results
mv ./chr${SGE_TASK_ID}_tmp_QC/tmp_passed_biallelic.txt tmp/chr${SGE_TASK_ID}_filtered_variants.txt
rm -r ./chr${SGE_TASK_ID}_tmp_QC

# Get genotype information at these SNPs and    
plink2 \
--extract tmp/chr${SGE_TASK_ID}_filtered_variants.txt \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${SGE_TASK_ID}_v3.bgen ref-first \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--threads 3 \
--memory 15000 \
--make-pgen \
--out tmp/chr${SGE_TASK_ID}_all_variants

# Create LD-pruned list of independent variants (+/- 500kb, steps of 50, R2 < 0.5)
plink2 \
--pfile tmp/chr${SGE_TASK_ID}_all_variants \
--indep-pairwise 500 kb 0.5 \
--threads 3 \
--memory 15000 \
--out chr${SGE_TASK_ID}

# Make pgen files for obesity-associated variants on these chromosomes only
plink2 \
--extract chr${SGE_TASK_ID}.prune.in \
--pfile tmp/chr${SGE_TASK_ID}_all_variants \
--threads 3 \
--memory 15000 \
--make-pgen \
--out pgen_format/chr${SGE_TASK_ID}_pruned

# Make bed files (hard-call thresholds needed for some functions)
plink2 \
--pfile pgen_format/chr${SGE_TASK_ID}_pruned \
--threads 3 \
--memory 15000 \
--make-bed \
--out bed_format/chr${SGE_TASK_ID}_pruned

# Remove genotype files for all variants (only retain pruned)
rm tmp/chr${SGE_TASK_ID}_all_variants.p*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
