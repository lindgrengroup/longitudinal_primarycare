#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 11/01/22

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1 
#$ -N ldsc_snp_heritability
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

cat /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt | while read STRATA || [[ -n $STRATA ]];
do
	# Get number of individuals in strata
	NSTRATA=$(wc -l < /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA}_ids.incl)
	
	# Quantitative traits from BOLT output
	for PARAMETER in "lmm_intercepts" "lmm_slopes_adj_baseline" "cspline_intercepts"
	do
		# Move to strata directory
		cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/${STRATA}/${PARAMETER}
		rm -r ldsc_h2
		mkdir ldsc_h2
		cd ldsc_h2

		# Check sumstats are in the right format
		munge_sumstats.py \
		--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_results/${STRATA}/${STRATA}_${PARAMETER}_final.txt.gz \
		--N ${NSTRATA} \
		--snp SNP --a1 Tested_Allele --a2 Other_Allele --p PVALUE \
		--frq AF_Tested --signed-sumstats BETA,0 \
		--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
		--out tmp

		# Run LDSC to get heritability 
		ldsc.py \
		--h2 tmp.sumstats.gz \
		--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--out ${STRATA}_${PARAMETER}_h2_ldsc

		rm tmp*
	done

	# Cluster traits from SAIGE output
	for KI in {1..6}
	do
		cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/post_GWAS/${STRATA}/k${KI}
		rm -r ldsc_h2
		mkdir ldsc_h2
		cd ldsc_h2

		# Check sumstats are in the right format
		munge_sumstats.py \
		--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_results/${STRATA}/${STRATA}_k${KI}_final.txt.gz \
		--N ${NSTRATA} \
		--snp SNP --a1 Tested_Allele --a2 Other_Allele --p PVALUE \
		--frq AF_Tested --signed-sumstats BETA,0 \
		--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
		--out tmp

		# Run LDSC to get heritability 
		ldsc.py \
		--h2 tmp.sumstats.gz \
		--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
		--out ${STRATA}_k${KI}_h2_ldsc

		rm tmp*
	done
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



