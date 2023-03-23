#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J extract_gws_snps
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS/logs/extract_gws_sig_snps-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	# Get list of all possible GWS SNPs from both b0 and average_trait GWAS
	awk -F '\t' 'NR>1 { if ($9 <= 5E-8) print $1 }' BOLT_results/${STRATA}_b0_final.txt \
	> post_GWAS/gws_sig_snps/tmp_snp_list_${STRATA}.txt
	awk -F '\t' 'NR>1 { if ($9 <= 5E-8) print $1 }' BOLT_results/${STRATA}_average_trait_final.txt \
	>> post_GWAS/gws_sig_snps/tmp_snp_list_${STRATA}.txt

	# Extract all results for these SNPs from each GWAS
	head -n 1 BOLT_results/${STRATA}_b0_final.txt > post_GWAS/gws_sig_snps/${STRATA}_b0_gws_sig_hits.txt
	grep -F post_GWAS/gws_sig_snps/tmp_snp_list_${STRATA}.txt BOLT_results/${STRATA}_b0_final.txt \
	>> post_GWAS/gws_sig_snps/${STRATA}_b0_gws_sig_hits.txt

	head -n 1 BOLT_results/${STRATA}_average_trait_final.txt > post_GWAS/gws_sig_snps/${STRATA}_average_trait_gws_sig_hits.txt
	grep -F post_GWAS/gws_sig_snps/tmp_snp_list_${STRATA}.txt BOLT_results/${STRATA}_average_trait_final.txt \
	>> post_GWAS/gws_sig_snps/${STRATA}_average_trait_gws_sig_hits.txt
done

rm post_GWAS/gws_sig_snps/tmp_snp_list_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0

