#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N prep_classification
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS

STRATA="BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb"

for STRATA_NAME in $STRATA
do
	mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants
	mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/reported_variants
	mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/condnl_results
	mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_ukbb_bed
	mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_ukbb_ld

	# TO RUN FOR ALL GENOME-WIDE SIGNIFICANT SNPS
	# Get significant (P <= 5E-8) SNPs 
	zcat BOLT_results/${STRATA_NAME}_lmm_intercepts_final.txt.gz | \
	awk -F '[[:space:]]+' 'NR>1 { if ($9 <= 5E-8) print $1,$2,$3 }' - \
	> post_GWAS/${STRATA_NAME}/lmm_intercepts_sig_snps.txt

	# Create temporary sumstats file in GCTA-COJO format and add column for sample size
	zcat BOLT_results/${STRATA_NAME}_lmm_intercepts_final.txt.gz | \
	awk -v n_gp="$N_GP" 'BEGIN{FS=OFS="\t"} {print $1, $4, $5, $6, $7, $8, $9, n_gp}' - \
	> post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_sumstats_gcta.txt
done 


echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
