#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 10/02/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J extract_gws_snps
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS/logs/extract_gws_sig_snps-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	for CLUSTER in k1 k1_k2 k1_k2_k3 k1_no_FU_adjustment k1_k2_no_FU_adjustment k1_k2_k3_no_FU_adjustment; do
		head -n 1 BOLT_results/${STRATA}_${CLUSTER}_final.txt > post_GWAS/lead_snps/all_strata_lead_snp_assocns/${STRATA}_${CLUSTER}_lead_snp_assocns.txt
		grep -F -w -f post_GWAS/lead_snps/all_strata_lead_snps.txt BOLT_results/${STRATA}_${CLUSTER}_final.txt \
		>> post_GWAS/lead_snps/all_strata_lead_snp_assocns/${STRATA}_${CLUSTER}_lead_snp_assocns.txt
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
