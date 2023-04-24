#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 12/04/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J lookup_snps_interest
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/lookup_snps_interest-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity
mkdir enrichment_change_signal_baseline

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	# Get list of SNPs
	awk 'NR > 1 {print $1}' 2211_models/GWAS/post_GWAS/lead_snps/${STRATA}_b0_final.lead_snps.txt \
	> enrichment_change_signal_baseline/${STRATA}_b0_lead_snps_rsids.txt
	# Get b1 sumstats
	head -n 1 2211_models/GWAS/BOLT_results/${STRATA}_b1_assoc.stats \
	> enrichment_change_signal_baseline/${STRATA}_b1_grepped_sumstats.txt
	grep -F -w -f enrichment_change_signal_baseline/${STRATA}_b0_lead_snps_rsids.txt \
	2211_models/GWAS/BOLT_results/${STRATA}_b1_assoc.stats \
	>> enrichment_change_signal_baseline/${STRATA}_b1_grepped_sumstats.txt
	# Get sumstats from clusters
	for CLUSTK in k1 k1_k2 k1_k2_k3; do
		head -n 1 highdim_splines/standardised_outcomes/GWAS/BOLT_results/${STRATA}_${CLUSTK}_assoc.stats \
		> enrichment_change_signal_baseline/${STRATA}_${CLUSTK}_grepped_sumstats.txt
		grep -F -w -f enrichment_change_signal_baseline/${STRATA}_b0_lead_snps_rsids.txt \
		highdim_splines/standardised_outcomes/GWAS/BOLT_results/${STRATA}_${CLUSTK}_assoc.stats \
		>> enrichment_change_signal_baseline/${STRATA}_${CLUSTK}_grepped_sumstats.txt
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
