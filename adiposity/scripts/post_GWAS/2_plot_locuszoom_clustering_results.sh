#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 19/04/21

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N weight_sexcomb_locuszoom_plots
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

for KI in k1 k1_k2_k3 k5 
do
	# Before running this, create batch list of SNPs to plot
	cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS

	# Create METAL-style input file for locuszoom
	zcat SAIGE_results/Weight_sex_comb/Weight_sex_comb_${KI}_final.txt.gz | \
	awk 'BEGIN{OFS="\t"} {print $1,$9}' - \
	> post_GWAS/Weight_sex_comb/${KI}/locuszoom/${KI}_input_sumstats.txt

	cd post_GWAS/Weight_sex_comb/${KI}/locuszoom

	# Run locuszoom
	module purge
	module load Python/2.7.18-GCCcore-10.2.0

	/well/lindgren/resources/locuszoom/bin/locuszoom \
	--metal ${KI}_input_sumstats.txt \
	--markercol SNP --pvalcol PVALUE \
	--hitspec /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/post_GWAS/Weight_sex_comb/lead_snps_input.txt \
	--pop EUR --build hg19 --source 1000G_March2012 
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
