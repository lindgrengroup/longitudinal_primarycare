#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 27/03/23

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 5
#SBATCH --array 1-22
#SBATCH -J get_longevity_rsids
#SBATCH -o /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/logs/make_bedfiles-%j.out

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity

grep "${SLURM_ARRAY_TASK_ID}_" gwascat_indep_longevity_rsids.txt | awk '{print $2}' > tmp_chr${SLURM_ARRAY_TASK_ID}_longevity_rsids.txt

plink2 \
--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v3.bgen \
--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
--extract tmp_chr${SLURM_ARRAY_TASK_ID}_longevity_rsids.txt \
--make-bed \
--out bedfiles/chr${SLURM_ARRAY_TASK_ID}_indep_rsids \
--memory 15000 \
--threads 5 

rm tmp_chr${SLURM_ARRAY_TASK_ID}_longevity_rsids.txt

# Once these have all run, concatenate with the following code that removes multi-allelic variants

module purge

cd bedfiles
touch tmp_files_to_concatenate.txt
for CHR in {1..22}; do
	echo "chr${CHR}_indep_rsids" >> tmp_files_to_concatenate.txt
done

plink2 \
--pmerge-list tmp_files_to_concatenate.txt bfile \
--pmerge-list-dir /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/sensitivity_longevity/bedfiles \
--max-alleles 2 \
--threads 1 \
--memory 15000 \
--make-bed \
--out longevity_indep_rsids

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
