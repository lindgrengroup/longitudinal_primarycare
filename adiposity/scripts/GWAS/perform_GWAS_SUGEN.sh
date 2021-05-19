#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 4
#$ -t 1-24 -tc 2
#$ -N adiposity_slope_SUGEN_GWAS
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Load SUGEN module
module load SUGEN/8.1-GCCcore-6.4.0

##########################################################################################

strata_name=`sed -n -e "$SGE_TASK_ID p" strata_filenames.txt`

# Make a new directory for GWAS results in the relevent strata
mkdir ${strata_name}

# Run SUGEN on each chromosome
for CHR_ID in {1..22}
do
	SUGEN \
	--pheno slope_files/${strata_name}.txt \
	--formula "RINTed_resids=genotyping_array" \
	--id-col IID \
	--family-col FID \
	--unweighted \
	--dosage \
	--vcf /well/lindgren/UKBIOBANK/samvida/general_resources/UKB_IMPUTED_VCF/ukb_imp_chr${CHR_ID}_v3.vcf.gz \
	--model linear \
	--hetero-variance FUyr_quartile \
	--out-prefix ${strata_name}/chr${CHR_ID}_assoc
done