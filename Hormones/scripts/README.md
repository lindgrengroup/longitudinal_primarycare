Scripts in this folder:

## GWAS
1. **1_prep_cross_sec.R** - Gather primary care hormone measures from pre-QCd data, see README in *../../full_primary_care/* for details on how this data was generated. Get cross-sectional value and age as earliest observed value closest to each individual's median. Split the data by sex and in each sex strata as well as sex-combined, adjust for covariates (age, age-squared), and rank-based inverse normal transform the values. 
2. **2_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre. Write QCd trait files. 
3. **3_perform_GWAS_BOLT.sh** - Array job (for all traits and strata) to perform GWAS under the linear mixed model framework in BOLT, using imputed genotypes from UK Biobank.
4. **4_filter_GWAS_results.R** - Filter GWAS results based on MAF (>0.1%), HWE pvalue (>1E-06), missingness (<5%), remove implausible standard errors (>10), and duplicate SNPs. Plot QQ-plots in different MAF bins and plot overall Manhattan plot for results.

## Phenotype ascertainment
1. **0_generate_controls.R** - Create list of potential controls for hormone ascertainment tests, which consists of individuals of White British ancestry with at least one record in GP data, but none of the 9 hormones tested.
2. **1a_phenome_wide_enrichment.R** - Fisher's test for enrichment of 307 binary disease codes (see *../../samvida_general/UKB_scripts/build_eid_phenotype_matrix.R* for details on how these were generated) in "cases" with a hormone measured in GP data vs "controls" generated above. Both sex-specific and sex-combined analyses performed. Important: this is time-agnostic, and disease diagnosis may be recorded at any time relative to hormone measurement.
3. **2a_plot_phenome_wide_enrichment.R** - Manhattan-style PheWAS plots for enrichment test p-values.

## Archived
1. **get_basic_stats.R** - Gathers primary care hormone measures; then plots basic information such as age distribution, sex distribution, etc.
2. **plot_trajectories.R** - Plots trajectories for individuals stratified by number of observations, overall change in hormone levels, etc. ADD MORE AS MORE WAYS TO LOOK AT DATA
3. **define_case_controls.R** - Defines case IDs and three strategies for control IDs: random population cohort, rando non-case cohort, and matched (age, BMI) with "cases" as individuals with specified hormone measure