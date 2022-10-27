#!/bin/bash

#SBATCH -A lindgren.prj
#SBATCH -p short
#SBATCH -c 1
#SBATCH -J get_signif_hits

echo "########################################################"
echo "Slurm Job ID: $SLURM_JOB_ID" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "##########################################################"

cd /well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity

rm ukb_no_gp/tmp.txt
touch ukb_no_gp/tmp.txt
for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	zcat 2204_models/GWAS/BOLT_results/${STRATA}_lmm_slopes_adj_int_final.txt.gz | \
	awk 'BEGIN {FS = "[[:space:]]+"} { if ($9<=5E-08) print $1,$2,$3 }' - \
	>> ukb_no_gp/tmp.txt
done
cat ukb_no_gp/tmp.txt | sort | uniq > ukb_no_gp/lmm_slopes_snps_to_replicate.txt

rm ukb_no_gp/tmp.txt
touch ukb_no_gp/tmp.txt
for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	for clustk in k1 k2 k3 k4; do
		awk 'BEGIN {FS = "[[:space:]]+"} { if ($9<=5E-08) print $1,$2,$3 }' \
		highdim_splines/GWAS/BOLT_results/${STRATA}_${clustk}_final.txt \
		>> ukb_no_gp/tmp.txt
	done
done
cat ukb_no_gp/tmp.txt | sort | uniq > ukb_no_gp/clustprobs_snps_to_replicate.txt

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	zcat 2204_models/GWAS/BOLT_results/${STRATA}_lmm_slopes_adj_int_final.txt.gz | \
	head -n 1 > 2204_models/GWAS/BOLT_results/${STRATA}_lmm_slopes_adj_int_gws_sig_hits.txt
	zcat 2204_models/GWAS/BOLT_results/${STRATA}_lmm_slopes_adj_int_final.txt.gz | \
	awk 'BEGIN {FS = "[[:space:]]+"} { if ($9<=5E-08 && $6>=0.01 && $6<=0.99) print $0 }' - \
	>> 2204_models/GWAS/BOLT_results/${STRATA}_lmm_slopes_adj_int_gws_sig_hits.txt
done

for STRATA in BMI_F BMI_M BMI_sex_comb Weight_F Weight_M Weight_sex_comb; do
	for clustk in k1 k2 k3 k4; do
		head -1 highdim_splines/GWAS/BOLT_results/${STRATA}_${clustk}_final.txt \
		> highdim_splines/GWAS/BOLT_results/${STRATA}_${clustk}_gws_sig_hits.txt
		awk '{ if ($9<=5E-08 && $6>=0.01 && $6<=0.99) print $0 }' \
		highdim_splines/GWAS/BOLT_results/${STRATA}_${clustk}_final.txt \
		>> highdim_splines/GWAS/BOLT_results/${STRATA}_${clustk}_gws_sig_hits.txt
	done
done

echo "###########################################################"
echo "Finished at: "`date` 
echo "###########################################################"
exit 0
