#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 2
#SBATCH -J replicate_associations

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

FILEPATH="" # REDACTED 

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	mkdir ${STRATA}
	mkdir tmp_${STRATA}

	# Run through each relevant chromosome
	# Need to specify VIF and max correlation because baseline age and age-sq are correlated covariates
	for CHR in {1..22}; do
		if [ -s ../replication_genotypes/chr${CHR}_snps.txt ]; then
			if [[ "${STRATA}" == "sex_comb" ]]; then 
				# Include sex flag as covariate

				# For lmm slopes
				plink2 \
				--pfile ${FILEPATH}/chr${CHR} \
				--pheno ${FILEPATH}/traits_for_GWAS/${STRATA}.txt \
				--pheno-name b1 \
				--covar-name genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--sex --covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/b1_chr${CHR}

				# For cluster probabilities
				plink2 \
				--pfile ${FILEPATH}/chr${CHR} \
				--pheno ${FILEPATH}/traits_for_GWAS/${STRATA}.txt \
				--pheno-name k1,k1_k2,k1_k2_k3 \
				--covar-name baseline_trait,baseline_age,age_sq,FU_n,FUyrs,year_of_birth,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--sex --covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/softkprobs_chr${CHR}
			else
				# Do not include sex flag as covariate
				plink2 \
				--pfile ${FILEPATH}/chr${CHR} \
				--pheno ${FILEPATH}/traits_for_GWAS/${STRATA}.txt \
				--pheno-name b1 \
				--covar-name genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/b1_chr${CHR}

				# For cluster probabilities
				plink2 \
				--pfile ${FILEPATH}/chr${CHR} \
				--pheno ${FILEPATH}/traits_for_GWAS/${STRATA}.txt \
				--pheno-name k1,k1_k2,k1_k2_k3 \
				--covar-name baseline_trait,baseline_age,age_sq,FU_n,FUyrs,year_of_birth,UKB_assmt_centre,genotyping.array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21 \
				--covar-variance-standardize --vif 1000 --max-corr 1 \
				--threads 3 --memory 15000 \
				--glm hide-covar \
				--out tmp_${STRATA}/softkprobs_chr${CHR}
			fi
		fi
	done 

	# Collate results from all chromosomes across traits
	head -n 1 tmp_${STRATA}/b1_chr1.b1.glm.linear > ${STRATA}/${STRATA}_b1.txt
	tail -n +2 -q tmp_${STRATA}/b1_chr*.b1.glm.linear >> ${STRATA}/${STRATA}_b1.txt

	for TRAIT in k1 k1_k2 k1_k2_k3; do
		head -n 1 tmp_${STRATA}/softkprobs_chr1.${TRAIT}.glm.linear > ${STRATA}/${STRATA}_${TRAIT}.txt
		tail -n +2 -q tmp_${STRATA}/softkprobs_*.${TRAIT}.glm.linear >> ${STRATA}/${STRATA}_${TRAIT}.txt
	done

	rm -r tmp_${STRATA}
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
