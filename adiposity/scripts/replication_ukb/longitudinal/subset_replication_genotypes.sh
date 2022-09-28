#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 4
#SBATCH -J extract_replication_genotypes

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

# This extracts SNPs to replicate
# But for ALL INDIVIDUALS
# When running GWAS, subset to the list of individuals we have phenotype data for

module load PLINK/2.00a2.3_x86_64

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp

for CHR in {1..22}; do
	rm tmp_varlist.txt
	touch tmp_varlist.txt
	awk -v CHR="$CHR" '{ if ($2==CHR) print $1 }' \
	lmm_slopes_snps_to_replicate.txt >> tmp_varlist.txt

	awk -v CHR="$CHR" '{ if ($2==CHR) print $1 }' \
	clustprobs_snps_to_replicate.txt >> tmp_varlist.txt

	cat tmp_varlist.txt | sort | uniq > replication_genotypes/chr${CHR}_snps.txt
	rm tmp_varlist.txt

	# Get genotype information at these SNPs
	if [ -s replication_genotypes/chr${CHR}_snps.txt ]; then
        plink2 \
		--extract replication_genotypes/chr${CHR}_snps.txt \
		--bgen /well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/ukb_imp_chr${CHR}_v3.bgen ref-first \
		--sample /well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample \
		--threads 3 \
		--memory 15000 \
		--make-pgen \
		--out replication_genotypes/chr${CHR}

		# Make bed files (hard-call thresholds needed for some functions)
		plink2 \
		--pfile replication_genotypes/chr${CHR} \
		--threads 3 \
		--memory 15000 \
		--make-bed \
		--out replication_genotypes/chr${CHR}
	fi
done

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
