cat strata_filenames.txt | while read STRATA_NAME || [[ -n $STRATA_NAME ]];
do
	# Get sample IDs to retain 
	awk 'NR > 1 {print $1}' /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA_NAME}_ids_passed_qc.txt > /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA_NAME}_ids_for_finemap.txt
	# Print number of samples 
	wc -l /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/sample_qc/${STRATA_NAME}_ids_for_finemap.txt
done
