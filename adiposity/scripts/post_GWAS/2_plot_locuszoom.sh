#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -t 1-7 -tc 2 
#$ -N adiposity_locuszoom_plots
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Before running this, create batch list of SNPs to plot

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" locuszoom_filenames.txt`

cd /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS

# Create METAL-style input file for locuszoom
zcat BOLT_filtered/${STRATA_NAME}_lmm_slopes_adj_baseline.txt.gz | \
awk 'BEGIN{OFS="\t"} {print $1,$9}' - \
> post_GWAS/${STRATA_NAME}/locuszoom/lmm_slopes_input_sumstats.txt

cd post_GWAS/${STRATA_NAME}/locuszoom

# Run locuszoom
module purge
module load Python/2.7.18-GCCcore-10.2.0

/well/lindgren/resources/locuszoom/bin/locuszoom \
--metal lmm_slopes_input_sumstats.txt \
--markercol SNP --pvalcol P_BOLT_LMM_INF \
--hitspec lead_snps_input.txt \
--pop EUR --build hg19 --source 1000G_March2012 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0