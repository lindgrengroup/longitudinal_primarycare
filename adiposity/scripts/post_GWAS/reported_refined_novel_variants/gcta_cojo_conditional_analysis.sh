#!/bin/bash

# Author: Samvida S. Venkatesh
# Date: 17/04/22

#$ -cwd
#$ -P lindgren.prjc -q short.qc
#$ -pe shmem 1
#$ -N gcta_cojo
#$ -j y

echo `date`: Executing task ${SGE_TASK_ID} of job ${JOB_ID} on `hostname` as user ${USER}

##########################################################################################

echo "passing covariates..."
covars=${covars//|/,}
echo ${covars}

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS/${STRATA}/classify_b0_variants

# Load PLINK module
module load PLINK/2.00a2.3_x86_64

# Read through all variants in file
while read var; do

	### PREPARE DATA

	# Create a temporary file for single variant	
	echo ${var} > tmp_${var}.txt
	# Remove SNP of interest from reported SNP list
	grep -vwE ${var} \
	reported_variants/${var}_published.txt > tmp_pub_around_${var}.txt	

	# Create temporary bed/bim/fam files for the relevant genomic region
	# Only retaining individuals of WB-ancestry for LD calculation
	plink2 \
	--bfile /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_IMPUTED_BED/ukb_imp_chr${CHR}_v3 \
	--keep /well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/fid_iid_white_british.txt \
	--snp ${var} \
	--window 1000 \
	--rm-dup force-first \
	--make-bed \
	--threads 1 \
	--memory 15000 \
	--out tmp_ukbb_bed/${var}_region

	# Create temporary LD-matrix for SNP vs published/reported SNPs in the region
	/apps/well/plink/1.90b3/plink \
	--bfile tmp_ukbb_bed/${var}_region \
	--r2 inter-chr \
	--ld-snp ${var} \
	--threads 1 \
	--memory 15000 \
	--out tmp_ukbb_ld/tmp_${var}_region
	# Only retain SNPs of interest (published)
	grep -f tmp_pub_around_${var}.txt \
	tmp_ukbb_ld/tmp_${var}_region.ld > tmp_ukbb_ld/${var}_region.ld
	rm tmp_ukbb_ld/tmp_${var}_region*

	### CONDITIONAL ANALYSIS

	# First remove collinear SNPs from conditional SNP-list
	/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
	--bfile tmp_ukbb_bed/${var}_region  \
	--cojo-file tmp_sumstats_gcta.txt \
	--extract reported_variants/${var}_published.txt \
	--cojo-slct \
	--out reported_variants/${var}_indpt_reported_snps

	# If there are no independent SNPs then the SNP of interest
	# is previously reported
	if [ ! -e reported_variants/${var}_indpt_reported_snps.jma.cojo ]; then
		echo ${var} >> reported_snp_list.txt
	else
		# Extract independent SNP list
		awk -F '\t' '(NR>1) {print $2}' \
		reported_variants/${var}_indpt_reported_snps.jma.cojo \
		> reported_variants/${var}_indpt_reported_snps.txt

		# If our SNP of interest is in the independently associated 
		# list of SNPs, then it may be novel or refined;
		# else it is reported
		if grep -q ${var} reported_variants/${var}_indpt_reported_snps.txt; then 
			# Remove SNP of interest from reported SNP list
			grep -vwE ${var} \
			reported_variants/${var}_indpt_reported_snps.txt \
			 > reported_variants/${var}_indpt_reported_snps_for_cond.txt	
		
			# If there are other SNPs in the independently associated list,
			# our SNP may be refined
			# else it is potentially novel
			if [ -s reported_variants/${var}_indpt_reported_snps_for_cond.txt ]; then
				echo ${var} >> potential_refined_snps.txt

				# 1. Effect of SNP of interest, conditioned on remaining SNPs
	  			/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
				--bfile tmp_ukbb_bed/${var}_region  \
				--cojo-file tmp_sumstats_gcta.txt \
				--extract reported_variants/${var}_indpt_reported_snps.txt \
				--cojo-cond reported_variants/${var}_indpt_reported_snps_for_cond.txt \
				--out condnl_results/${var}_cond_on_pub
				
				# 2. Effect of published SNPs, conditioned on variant of interest
	  			/apps/well/gcta/1.91.5beta/gcta_1.91.5beta/gcta64 \
				--bfile tmp_ukbb_bed/${var}_region  \
				--cojo-file tmp_sumstats_gcta.txt \
				--extract reported_variants/${var}_published.txt \
				--cojo-cond tmp_${var}.txt \
				--out condnl_results/pub_cond_on_${var}

			else
				echo ${var} >> potential_novel_snps.txt
			fi
		else
			echo ${var} >> reported_snp_list.txt
		fi
	fi
	rm tmp_${var}.txt
	rm tmp_pub_around_${var}.txt
	rm tmp_ukbb_bed/${var}_region*
	rm reported_variants/${var}_indpt_*
done < tmp_variant_list_chr${CHR}.txt 

echo "###########################################################"
echo "Finished at: "`date`
echo "###########################################################"
exit 0
