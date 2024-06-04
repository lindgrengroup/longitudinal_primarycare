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

module load Python/3.10.8-GCCcore-12.2.0
source /well/lindgren/samvida/python/gwas-catalog-sumstats-${MODULE_CPU_TYPE}/bin/activate

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity

for ADIPO in BMI Weight; do
	for SS in F M sex_comb; do
		# LMM GWASs
		for LMM in b0 b1; do
			# old file: 1.SNP 2.CHR 3.POS 4.Tested_Allele 5.Other_Allele 6.AF_Tested 7.BETA 8.SE 9.PVALUE
			# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id
			awk -F'\t' -v OFS='\t' '{print $2, $3, $4, $5, $7, $8, $6, $9, $1}' 2211_models/GWAS/BOLT_results/${ADIPO}_${SS}_${LMM}_final.txt \
			> sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_sumstats.tsv
			sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id/' \
			sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_sumstats.tsv

			# Validate
			gwas-ssf validate --errors-out sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_sumstats.tsv
			# gzip
			gzip sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${LMM}_sumstats.tsv
		done

		# Cluster belonging GWASs
		for CLUST in k1 k1_k2 k1_k2_k3; do
			# old file: 1.SNP 2.CHR 3.POS 4.Tested_Allele 5.Other_Allele 6.AF_Tested 7.BETA 8.SE 9.PVALUE
			# new file: chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rs_id
			awk -F'\t' -v OFS='\t' '{print $2, $3, $4, $5, $7, $8, $6, $9, $1}' highdim_splines/standardised_outcomes/GWAS/BOLT_results/${ADIPO}_${SS}_${CLUST}_final.txt \
			> sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_sumstats.tsv
			sed -i '1s/.*/chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\teffect_allele_frequency\tp_value\trs_id/' \
			sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_sumstats.tsv

			# Validate
			gwas-ssf validate --errors-out sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_sumstats.tsv
			# gzip
			gzip sumstats_gwas_catalog_deposit/${ADIPO}_${SS}_${CLUST}_sumstats.tsv
		done
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
