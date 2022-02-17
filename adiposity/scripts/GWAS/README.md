Scripts in this folder:

1. **1_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre. This only needs to be run once for each strata, as all models share the same individuals.
2. **2_RINT_traits.R** - Get GWAS phenotype, i.e. adjust the spline or LMM coefficients for covariates and rank-based inverse normal transform the residuals within each sex strata.\
**2a_RINT_traits_within_b0_quartiles.R** - Perform same adjustments and transformations as (2) within each quartile of intercepts from linear mixed models.\
**2b_collate_traits_PLINK.R** - Format phenotype files for PLINK software association testing.
**2c_create_cluster_GWAS_files.R** - Format phenotype files as cluster membership and covariates for association testing in SAIGE with binary cluster membership.
3. **3b_perform_variant_assocns_PLINK.sh** - Array job (for all traits and strata) to test for associations of slopes, intercepts, or other model parameters with a subset of metabolic and endocrine variants (see ../general_scripts/extract_gwascat_variants.R for how this was generated) under the generalised linear model framework in PLINK. \
**3c_perform_variant_assocns_b0_quartiles_PLINK.sh** - Same as (3b) but phenotypes are split into quartiles based on intercepts from linear mixed models, as found in (2a). 
4. **4b_ and 4c_plot_qq_manhattan.R** - QQ plots and Manhattan plots for associations of various model parameters with subset of metabolic and endocrine variants (as described in *../general_scripts/extract_gwascat_variants.R*) \
**generic_submit_BOLT_scripts.R** - Loops over all strata and parameters (lmm intercepts, lmm slopes adjusted for baseline, cubic spline intercepts) to submit BOLT GWAS, filtering, etc. Relies on helper scripts in *./BOLT_helpers/*. See README in *./BOLT_helpers/* for more details. \
**generic_submit_SAIGE_scripts.R** - Loops over all strata and clusters to submit SAIGE step 1, step 2, or results filtering for cluster membership as the GWAS trait. Relies on helper scripts in *./SAIGE_helpers/*. See README in *./SAIGE_helpers/* for more details.
