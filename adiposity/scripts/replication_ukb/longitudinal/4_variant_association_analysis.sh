#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J replicate_associations

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/longit_replication_GWAS
	mkdir ${STRATA}
	mkdir tmp_${STRATA}

	# Run through each relevant chromosome
	# Need to specify VIF and max correlation because baseline age and age-sq are correlated covariates
	for CHR in {1..22}; do
		if [ -s ../replication_genotypes/chr${CHR}_snps.txt ]; then
			if [[ "${STRATA}" == "sex_comb" ]]; then 
				# Include sex flag as covariate
				plink2 \
				--pfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/chr${CHR} \
				--pheno /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/traits_for_GWAS/${STRATA}.txt \
				--pheno-name lmm_slopes_adj_int,k1,k2,k3,k4 \
				--covar-name baseline_age,age_sq,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--sex --covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/chr${CHR}
			else
				# Do not include sex flag as covariate
				plink2 \
				--pfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/chr${CHR} \
				--pheno /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/replication_genotypes/traits_for_GWAS/${STRATA}.txt \
				--pheno-name lmm_slopes_adj_int,k1,k2,k3,k4 \
				--covar-name baseline_age,age_sq,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/chr${CHR}
			fi
		fi
	done 

	# Collate results from all chromosomes across traits

	for TRAIT in lmm_slopes_adj_int k1 k2 k3 k4; do
		head -n 1 tmp_${STRATA}/chr1.${TRAIT}.glm.linear > ${STRATA}/${STRATA}_${TRAIT}.txt
		tail -n +2 -q tmp_${STRATA}/chr*.${TRAIT}.glm.linear >> ${STRATA}/${STRATA}_${TRAIT}.txt
	done

	rm -r tmp_${STRATA}
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0