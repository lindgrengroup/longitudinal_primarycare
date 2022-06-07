#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 11/01/22

# Pan-UKBB LD panel (EUR) for 1 million HapMap variants

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -t 1-6 -tc 3 
#$ -N hrg_gp_vs_giant
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

STRATA_NAME=`sed -n -e "$SGE_TASK_ID p" bmi_whr_strata_filenames.txt`
GP_N=`sed -n -e "$SGE_TASK_ID p" gp_gwas_samplesizes.txt`
META_FILENAME=`sed -n -e "$SGE_TASK_ID p" giant_ukb_meta_filenames.txt`

cd /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/post_GWAS/${STRATA_NAME}
rm -r ldsc
mkdir ldsc
cd ldsc

module load LDSC/1.0.1-20200216-foss-2018b-Python-2.7.15

# Check sumstats are in the right format
munge_sumstats.py \
--sumstats /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/BOLT_filtered/${STRATA_NAME}_lmm_intercepts.txt.gz \
--N ${GP_N} \
--snp SNP --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM_INF \
--frq A1FREQ --signed-sumstats BETA,0 \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out tmp_gp

munge_sumstats.py \
--sumstats ${META_FILENAME} \
--snp SNP --N-col N --a1 Tested_Allele --a2 Other_Allele \
--p P --frq Freq_Tested_Allele --signed-sumstats BETA,0 --info INFO \
--merge-alleles /well/lindgren/resources/sldsc/aux_files/w_hm3.snplist \
--out tmp_meta

# Run LDSC to get heritability and genetic correlation
ldsc.py \
--rg tmp_gp.sumstats.gz,tmp_meta.sumstats.gz \
--ref-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--w-ld /well/lindgren/resources/panukb_ld_scores_complete/ldsc_format/rsids/UKBB.EUR.rsid \
--out h2_rg_giant

rm tmp_*

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0



