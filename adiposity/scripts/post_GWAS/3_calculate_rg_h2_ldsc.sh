#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 06/12/22

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 2
#$ -N hrg_b0_vs_change
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

# Need to pass: STRATA, SAMPLE_SIZE
echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/rg_h2
mkdir ${STRATA}
cd ${STRATA}

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# Check sumstats are in the right format
# b0 sumstats
munge_sumstats.py \
--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/BOLT_results/${STRATA}_b0_final.txt \
--N ${SAMPLE_SIZE} \
--snp SNP --a1 Tested_Allele --a2 Other_Allele --p PVALUE \
--frq AF_Tested --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out tmp_b0

# b1 sumstats
munge_sumstats.py \
--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/BOLT_results/${STRATA}_b1_final.txt \
--N ${SAMPLE_SIZE} \
--snp SNP --a1 Tested_Allele --a2 Other_Allele --p PVALUE \
--frq AF_Tested --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out tmp_b1

# Run LDSC to get heritability and genetic correlation
ldsc.py \
--rg tmp_b0.sumstats.gz,tmp_b1.sumstats.gz \
--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--out ${STRATA}_h2_rg_b0_b1

for clustk in k1 k1_k2 k1_k2_k3; do
	munge_sumstats.py \
	--sumstats /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/BOLT_results/${STRATA}_${clustk}_final.txt \
	--N ${SAMPLE_SIZE} \
	--snp SNP --a1 Tested_Allele --a2 Other_Allele --p PVALUE \
	--frq AF_Tested --signed-sumstats BETA,0 \
	--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
	--out tmp_${clustk}

	# Run LDSC to get heritability and genetic correlation
	ldsc.py \
	--rg tmp_b0.sumstats.gz,tmp_${clustk}.sumstats.gz \
	--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
	--out ${STRATA}_h2_rg_b0_${clustk}
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
