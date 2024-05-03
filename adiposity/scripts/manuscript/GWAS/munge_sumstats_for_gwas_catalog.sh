#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 02/05/24

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J munge_sumstats_gwascat
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sumstats_gwas_catalog_deposit/logs/munge-sumstats-%j.out

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity

for ADIPO in BMI Weight; do
	for SS in F M sex_comb; do
		# LMM GWASs
		for LMM in b0 b1; do
			# old file: 1.SNP 2.CHR 3.BP 4.GENPOS 5.ALLELE1 6.ALLELE0 7.A1FREQ 8.F_MISS 9.CHISQ_LINREG 10.P_LINREG 11.BETA 12.SE 13.CHISQ_BOLT_LMM_INF 14.P_BOLT_LMM_INF
			# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id
			zcat 2211_models/GWAS/BOLT_results/${ADIPO}_${SS}_${LMM}_filtered.txt.gz \
			| awk -F'\t' -v OFS='\t' '{print $2, $3, $5, $6, $11, $12, $7, $14, $1}' \
			> sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_filtered.tsv
			sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id/' \
			sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_filtered.tsv
			gzip sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_filtered.tsv
		done

		# Cluster belonging GWASs
		for CLUST in k1 k1_k2 k1_k2_k3; do
			# old file: 1.SNP 2.CHR 3.BP 4.GENPOS 5.ALLELE1 6.ALLELE0 7.A1FREQ 8.F_MISS 9.CHISQ_LINREG 10.P_LINREG 11.BETA 12.SE 13.CHISQ_BOLT_LMM_INF 14.P_BOLT_LMM_INF
			# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id
			zcat highdim_splines/standardised_outcomes/GWAS/BOLT_results/${ADIPO}_${SS}_${CLUST}_filtered.txt.gz \
			| awk -F'\t' -v OFS='\t' '{print $2, $3, $5, $6, $11, $12, $7, $14, $1}' \
			> sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_filtered.tsv
			sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id/' \
			sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_filtered.tsv
			gzip sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_filtered.tsv
		done
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
