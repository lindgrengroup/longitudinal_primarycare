#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J grep_longevity_rsids
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/logs/grep_longevity_rsids-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity

# Longevity rsids to grep
awk '{ print $2 }' sensitivity_longevity/gwascat_indep_longevity_rsids.txt > sensitivity_longevity/indep_rsids_only.txt

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	head -n 1 2211_models/GWAS/BOLT_results/${STRATA}_b1_final.txt > sensitivity_longevity/grepped_sumstats/${STRATA}_b1_longevity_rsids.txt
	grep -F -w -f sensitivity_longevity/indep_rsids_only.txt 2211_models/GWAS/BOLT_results/${STRATA}_b1_final.txt \
	>> sensitivity_longevity/grepped_sumstats/${STRATA}_b1_longevity_rsids.txt

	for CLUSTK in k1 k1_k2 k1_k2_k3; do
		head -n 1 highdim_splines/standardised_outcomes/GWAS/BOLT_results/${STRATA}_${CLUSTK}_final.txt \
		> sensitivity_longevity/grepped_sumstats/${STRATA}_${CLUSTK}_longevity_rsids.txt
		grep -F -w -f sensitivity_longevity/indep_rsids_only.txt highdim_splines/standardised_outcomes/GWAS/BOLT_results/${STRATA}_${CLUSTK}_final.txt \
		>> sensitivity_longevity/grepped_sumstats/${STRATA}_${CLUSTK}_longevity_rsids.txt
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
