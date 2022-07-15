#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N create_list_cond_snps
#$ -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/
#$ -e /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS

mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants
mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/reported_variants
mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/condnl_results
mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_ukbb_bed
mkdir post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_ukbb_ld

# TO RUN FOR ALL GENOME-WIDE SIGNIFICANT SNPS
# Get significant (P <= 5E-8) SNPs 
zcat BOLT_results/${STRATA_NAME}_lmm_intercepts_final.txt.gz | \
awk -F '[[:space:]]+' 'NR>1 { if ($9 <= 5E-8) print $1,$2,$3 }' - \
> post_GWAS/${STRATA_NAME}/lmm_intercepts_sig_snps.txt

# Create temporary sumstats file in GCTA-COJO format and add column for sample size
zcat BOLT_results/${STRATA_NAME}_lmm_intercepts_final.txt.gz | \
awk -v n_gp="$N_GP" 'BEGIN{FS=OFS="\t"} {print $1, $4, $5, $6, $7, $8, $9, n_gp}' - \
> post_GWAS/${STRATA_NAME}/classify_lmm_intercept_variants/tmp_sumstats_gcta.txt

# Create list of reported SNPs within +/- 500kb of every unreported, significant SNP
# and run GCTA-COJO for conditional analysis by looping over all variants
Rscript /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/classify_reported_novel_variants.R \
--strata=${STRATA_NAME}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
