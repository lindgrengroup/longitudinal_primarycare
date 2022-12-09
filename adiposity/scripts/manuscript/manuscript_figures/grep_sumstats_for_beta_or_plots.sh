#!/bin/bash

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS

for PHENO in BMI Weight; do
	for SEX_STRATA in F M sex_comb; do
		for CLUSTK in k1 k1_k2 k1_k2_k3; do
			cat BOLT_results/${PHENO}_${SEX_STRATA}_${CLUSTK}_final.txt | head -n 1 > post_GWAS/lead_snps/${PHENO}_${SEX_STRATA}_${CLUSTK}_all_leadsnp_assocns.txt
			grep -wFf all_lead_snps.txt BOLT_results/${PHENO}_${SEX_STRATA}_${CLUSTK}_final.txt >> post_GWAS/lead_snps/${PHENO}_${SEX_STRATA}_${CLUSTK}_all_leadsnp_assocns.txt
		done
	done
done
