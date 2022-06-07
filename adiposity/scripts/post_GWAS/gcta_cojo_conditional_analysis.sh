#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N gcta_cojo
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants

# Load PLINK module
module load PLINK/2.00a2.3_x86_64

while read var; do
  # Create temporary bed/bim/fam files for the relevant genomic region
	plink2 \
	--bfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_IMPUTED_BED/ukb_imp_chr${CHR}_v3 \
	--snp ${var} \
	--window 1000 \
	--make-bed \
	--threads 1 \
	--memory 15000 \
	--out tmp_ukbb_bed/${var}_region

# Run GCTA-COJO for conditional analysis

# First remove collinear SNPs from conditional SNP-list
	/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
	--bfile tmp_ukbb_bed/${var}_region  \
	--cojo-file tmp_sumstats_gcta.txt \
	--extract reported_variants/${var}_published.txt \
	--cojo-slct \
	--out reported_variants/${var}_indpt_reported_snps

	# Extract independent SNP list
	awk -F '\t' '(NR>1) {print $2}' \
	reported_variants/${var}_indpt_reported_snps.jma.cojo \
	> reported_variants/${var}_indpt_reported_snps.txt

	# Check if the original SNP is still in the non-collinear list
	# If it is, perform conditional analysis
	if grep -q ${var} reported_variants/${var}_indpt_reported_snps.txt; then
	    /apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
		--bfile tmp_ukbb_bed/${var}_region  \
		--cojo-file tmp_sumstats_gcta.txt \
		--extract reported_variants/${var}_indpt_reported_snps.txt \
		--cojo-cond reported_variants/${var}_indpt_reported_snps.txt \
		--out condnl_results/${var}
	else
	    echo ${var} >> reported_snp_list.txt
	    rm reported_variants/${var}_indpt_reported_snps.*.cojo
	fi
	rm tmp_ukbb_bed/${var}_region.*
done < /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/tmp_variant_list_chr${CHR}.txt

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
