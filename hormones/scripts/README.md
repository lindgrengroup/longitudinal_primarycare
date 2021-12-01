Scripts in this folder:

## GWAS

1. **1_prep_cross_sec.R** - Get cross-sectional hormone data from GP records by either extracting the single available value or (in the case of multiple measurements) extracting the first reported value closest to the individual's median hormone measurement. In sex-specific and sex-combined strata, adjust for covariates and rank-based inverse normal transform the residuals as phenotype for GWAS. 
2. **2_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre. 
3. **3_perform_GWAS_BOLT.sh** - Array job (for all hormones and strata) to perform GWAS under the linear mixed model framework in BOLT, using imputed genotypes from UK Biobank. 
4. **4_filter_GWAS_results.R** - Filter GWAS results based on MAF (>0.1%), HWE pvalue (>1E-06), missingness (<5%), remove implausible standard errors (>10), and duplicate SNPs. Plot QQ-plots in different MAF bins and plot overall Manhattan plot for results.

## PHENOTYPE ASCERTAINMENT
1. **0_generate_controls.R** - Create set of control IDs of White British ancestry that have at least one GP record, but none of the 9 hormones studied measured. 
2. **1a_phenome_wide_enrichment.R** - For each hormone, perform Fisher's test for enrichment of 307 binary diseases (see *../../samvida_general/UKB_scripts/build_eid_phenotype_matrix.R* for details on how these were obtained). Cases are defined as individuals with hormone measured in GP, controls as above. Performed in both sex-specific and sex-combined strata.
3. **2a_plot_phenome_wide_enrichment.R** - Manhattan-style PheWAS plots to plot p-values for enrichment of diseases (grouped by ICD code chapter) across all hormones and in hormone-specific plots.

## ARCHIVED

Clustering and trajectory plotting scripts that might be useful later.
