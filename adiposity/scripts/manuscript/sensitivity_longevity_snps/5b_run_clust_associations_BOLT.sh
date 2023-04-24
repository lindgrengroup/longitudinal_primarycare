#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 10
#SBATCH -J clust_assocns_longevity
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/logs/clust_assocns_longevity-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

/apps/well/bolt-lmm/2.3.2/bolt \
--bfile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/bedfiles/longevity_indep_rsids \
--phenoFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/traits_for_GWAS/${STRATA}_soft_clust_probs.txt \
--phenoCol=${CLUST} \
--covarFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/standardised_outcomes/GWAS/sample_qc/${STRATA}_covariates.txt \
--qCovarCol=baseline_trait \
--qCovarCol=baseline_age \
--qCovarCol=age_sq \
--qCovarCol=FUyrs \
--qCovarCol=FU_n \
--qCovarCol=PC{1:21} \
--covarCol=genotyping.array \
--covarCol=UKB_assmt_centre \
--covarCol=sex \
--covarMaxLevels=30 \
--maxMissingPerSnp 0.05 \
--LDscoresFile=/apps/well/bolt-lmm/2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--LDscoresMatchBp \
--geneticMapFile=/apps/well/bolt-lmm/2.3.2/tables/genetic_map_hg19_withX.txt.gz \
--lmm \
--numThreads=8 \
--statsFile=/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/BOLT_results/${STRATA}_${CLUST}_assoc.stats.gz \
--verboseStats

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
