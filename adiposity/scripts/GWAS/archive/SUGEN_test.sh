#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 29/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 3
#$ -N SUGEN_phenotype_BMI
#$ -o /well/lindgren/UKBIOBANK/samvida/adiposity/GWAS
#$ -j y

echo "########################################################"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started date: "`date`
echo "##########################################################"

module load SUGEN/8.1-GCCcore-6.4.0
mkdir test_slopes_noQC

SUGEN \
	--pheno slope_files/test_slopes.txt \
	--formula "RINTed_resids=genotyping_array" \
	--id-col IID \
	--family-col FID \
	--unweighted \
	--dosage \
	--vcf /well/lindgren/UKBIOBANK/samvida/general_resources/UKB_IMPUTED_VCF/ukb_imp_chr21_v3.vcf.gz \
	--model linear \
	--hetero-variance FUyr_quartile \
	--out-prefix test_slopes_noQC/chr21_assoc

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 