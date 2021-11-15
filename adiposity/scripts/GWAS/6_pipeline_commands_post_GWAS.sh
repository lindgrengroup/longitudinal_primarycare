cd /well/lindgren/UKBIOBANK/samvida/general_resources/pipeline

source activate pipeline

cat /well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/strata_filenames.txt | while read STRATA_NAME || [[ -n $STRATA_NAME ]];
do
	pipeline.r --outname=${STRATA_NAME}_lmm_intercepts \
	--outdir=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/post_GWAS \
	--input-json=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/post_GWAS_pipeline_inputs/${STRATA_NAME}_inputs.json \
	--job-dependencies-json=./json/job_dependencies.json \
	--sumstats=/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/BOLT_filtered/${STRATA_NAME}_lmm_intercepts.txt
done

