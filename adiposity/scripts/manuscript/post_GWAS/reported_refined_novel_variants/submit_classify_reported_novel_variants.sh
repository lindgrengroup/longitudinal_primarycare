#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N classify_reported_novel_variants
#$ -o post_GWAS/
#$ -e post_GWAS/
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

mkdir post_GWAS/${STRATA_NAME}/classify_b0_variants
mkdir post_GWAS/${STRATA_NAME}/classify_b0_variants/reported_variants
mkdir post_GWAS/${STRATA_NAME}/classify_b0_variants/condnl_results
mkdir post_GWAS/${STRATA_NAME}/classify_b0_variants/tmp_ukbb_bed
mkdir post_GWAS/${STRATA_NAME}/classify_b0_variants/tmp_ukbb_ld

# TO RUN FOR ALL GENOME-WIDE SIGNIFICANT SNPS
# Get significant (P <= 5E-8) SNPs 
# zcat BOLT_results/${STRATA_NAME}_lmm_intercepts_final.txt.gz | \
# awk -F '[[:space:]]+' 'NR>1 { if ($9 <= 5E-8) print $1,$2,$3 }' - \
# > post_GWAS/${STRATA_NAME}/lmm_intercepts_sig_snps.txt

# Create temporary sumstats file in GCTA-COJO format and add column for sample size
awk -v n_gp="$N_GP" 'BEGIN{FS=OFS="\t"} {print $1, $4, $5, $6, $7, $8, $9, n_gp}' BOLT_results/${STRATA_NAME}_b0_final.txt \
> post_GWAS/${STRATA_NAME}/classify_b0_variants/tmp_sumstats_gcta.txt

# Create list of reported SNPs within +/- 500kb of every unreported, significant SNP
# and run GCTA-COJO for conditional analysis by looping over all variants
Rscript scripts/classify_reported_novel_variants.R \
--strata=${STRATA_NAME}

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
