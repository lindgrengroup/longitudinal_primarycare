Scripts in this folder:

## MAIN ANALYSES

### Linear mixed models
1. **1_apply_lmms.R** - Build linear mixed effects models for effect of time (from baseline measurement) on trait, allowing for fixed and random intercepts and slopes; adjust for covariates (baseline age, age-squared, sex, genetic PCs, and data provider). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for intercept and slope as the fixed + random effect for each individual. 

### Cubic splines
1. **1_apply_splines.R** - Build cubic spline mixed effects models for effect of age on trait, allowing for fixed and random intercepts and splines. Run through models with increasing number of knots in the fixed effects (3-10) and random effects (1 - # fixed effects) to choose model with lowest BIC. Adjust for covariates (sex, genetic PCs, and data provider). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for each term as the fixed + random effect for each individual. 

### Cross-sectional (baseline)
1. **6_baseline_clusters.R** - Cluster the adjusted baseline trait value for each individual, stepping through 1:6 clusters allowing for varying shape, variance, etc. Assign individuals to clusters.


## PLOTTING AND DESCRIPTIVE STATISTICS

### Raw data - 
1. **report_characteristics.R** - Generate tables and figures for various descriptive characteristics of raw data.
2. **disease_associated_trajectories.R** - For a sample of diseases, plot raw data trajectories in individuals with and without disease.
### Modelling results - 
Change these scripts for each set of models, linear, cubic splines, etc.) - 
1. **plot_BLUP_distributions.R** - Plot randomly selected subset of BLUPs from models to check distributions as well as relationships between BLUPs and model covariates.
2. **plot_model_predictions.R** - Plot individual-level model predictions for randomly selected individuals, those with ve3ry few or many repeat measures, at the tails of model BLUPs, etc.  

## GWAS

1. **1_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre. This only needs to be run once for each strata, as all models share the same individuals.
2. **2_RINT_traits.R** - Get GWAS phenotype, i.e. adjust the spline or LMM coefficients for covariates and rank-based inverse normal transform the residuals within each sex strata.\
**2a_RINT_traits_within_b0_quartiles.R** - Perform same adjustments and transformations as (2) within each quartile of intercepts from linear mixed models.\
**2b_collate_traits_PLINK.R** - Format phenotype files for PLINK software association testing.
3. **3_perform_GWAS_BOLT.sh** - Array job (for all traits and strata) to perform GWAS under the linear mixed model framework in BOLT, using imputed genotypes from UK Biobank. Perform initial filtering by MAF > 1%, HWE p-value > 1E-6, INFO > 0.8, missingness < 5%, and bi-allelic SNPs. \
**3b_perform_variant_assocns_PLINK.sh** - Array job (for all traits and strata) to test for associations of slopes, intercepts, or other model parameters with a subset of metabolic and endocrine variants (see ../general_scripts/extract_gwascat_variants.R for how this was generated) under the generalised linear model framework in PLINK. \
**3c_perform_variant_assocns_b0_quartiles_PLINK.sh** - Same as (3b) but phenotypes are split into quartiles based on intercepts from linear mixed models, as found in (2a). 
4. **4_filter_GWAS_results.R** - Filter GWAS results on implausible standard errors (>10) and duplicate SNPs. Plot QQ-plots in different MAF bins and plot overall Manhattan plot for results.\
**4b_ and 4c_plot_qq_manhattan.R** - QQ plots and Manhattan plots for associations of various model parameters with subset of metabolic and endocrine variants (as described in *../general_scripts/extract_gwascat_variants.R*)

## PHENOTYPE ENRICHMENT
1. **female_phenotype_enrichment.R** - Get age at menarche, age at menopause, and number of live births from UKB to compare distributions in different clusters in female-specific analyses.
2. **cluster_disease_associations.R** - Use logistic regression (with dummy variables for each cluster), adjusted for covariates, to test for association of clusters with >300 diseases identified by ICD and primary care codes. Calculate prevalence of each disease in the cohort and each cluster. Alternative: Fisher's enrichment test for disease counts in each cluster vs cohort.

## MULTIVARIATE ANALYSES, PLOTS, ENRICHMENT, ETC.
1. **focus_BMI_WHR_overlap.R** - Data summaries and enrichment tests for individuals with both BMI and WHR trajectories in the data.
2. **multivariate_report_characteristics.R** - Data summaries, overlaps, and trajectories for individuals who have more than one adiposity trait measured in the data.

## ARCHIVED:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Removes population-level noise by calculating thresholds based on UKBIOBANK-measured BMI and trajectories of extreme individuals 
3. **bmi_post_stage1_stats.R** - Plots basic information such as age distribution, sex distribution, longitudinal traits, etc. of cleaned BMI measures; same script used after individual outlier removal
4. **bmi_individual_outlier_removal.R** - Removes individual-level noise by calculating sequential logFC (fold-change between consecutive measurements) as well as accounting for time between consecutive measurements; individuals in the tails of these distributions are highlighted and the noisy observation (farthest from individual median or UKBIOBANK BMI) is removed
5. **bmi_trajectories.R** - Plots trajectories for individuals stratified by number of observations, overall change in BMI, etc.
6. **compare_adj_slope_models.R** - ANOVA, AIC, etc. to compare models for BLUP adjustment, plot mean trajectories by different model quartiles to compare effect of adjustment. 
7. **1_get_phenotypes.R** - Gather primary care trait measurements and add same traits from UKBIOBANK, removes individuals without repeat measures.
8. **2_clean_phenotypes.R** - For each trait, QC on age, implausible and extreme values, and remove individuals without repeat measures. Flags pregnancy. Inter-convert all BMI and weight measures.
9. **3_plot_QC.R** - Visual inspection for unrealistically large jumps between time-points. Remove measurements causing unrealistic jump.
10. **4_calculate_covariates.R** - Calculate individual-level phenotyping (baseline age, baseline BMI, follow-up years, etc.) and genotyping (array, PCs) covariates. 
11. **5_adjust.R** - Adjust traits for fixed effect covariates, such as baseline age, follow-up length, sex, genetic PCs, data provider, etc. Save adjusted and fitted values.

### Linear mixed models - ARCHIVED
1. **5_raw_slopes.R** - Perform mixed effects regression with fixed and random effects for individual and age, regressing adiposity trait on age. Calculate BLUP for raw slope for each individual, remove outliers > 5 S.D. away from mean in each stratum. 
2. **6_adjust_slopes.R** - Run linear models to adjust raw BLUPs for a range of covariates with nested models. Calculate variance explained by baseline model covariates and compare models with ANOVA. (6a) picks the best model to carry forward.
3. **7_adjust_slope_groups.R** - Calculate groups (gainers, quartiles, etc.).

### Mixture models - ARCHIVED
1. **5_baseline_models.R** - Run linear and spline mixed effects models in each sex- and ancestry-specific stratum (as above) with k = 1 to establish baseline model before increasing k. Outcome is adjusted trait value (adjusted for baseline age, baseline age-squared, baseline adiposity trait, number of follow-up years and number of follow-up measures, first 21 genetic principal components, and sex in sex-combined analyses) Calculate BIC to compare linear and spline models.
2. **6_stepwise_models.R** - After picking best baseline model (linear), run stepflexmix to go through k = 1:10. Plot BIC for each k-value and save the best model (with lowest BIC).

### Spline mixed models - ARCHIVED
1. **6_slopes.R** - Perform mixed effects regression for each trait on linear and non-linear (spline) terms for age in the fixed and random effects, with the best model calculated by BIC. Calculate random + fixed effect coefficients for age-related terms for each individual, remove outliers > 5 S.D. away from mean in each term. Plot model predictions. Parses trait as command line argument, so submit jobs with **6_submit_slopes.sh**
2. **7_find_clusters.R** - Cluster the fixed + random effect coefficients generated in the previous script with mixture model clustering, stepping through 1:6 clusters allowing for varying shape, variance, etc. Assign individuals to clusters.

### GWAS - ARCHIVED
1. **prepare_GWAS_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Get GWAS phenotype - RINTed adjusted BLUP - for QCd individuals in full cohort as well as in only gainers.
2. **perform_GWAS_{BOLT/SUGEN}.sh** - Perform genome-wide association study using specified software.

### Regression spline mixed models - ARCHIVED
5. **5a_test_increasing_degrees.R** and **5b_apply_regression_models.R** - Run polynomial spline effects of age (fixed and random effect) for increasing degrees (1:10) and plot heatmap of fixed effect coefficients to determine number of degrees in final model. Once the polynomial degree is chosen (cubic), run natural cubic spline regression, adjusted for covariates, and save the random effect coefficients in each strata.
6. **6_random_effect_distributions.R** - Plot coefficients of the random effect terms to check their distributions, and also test for association with baseline covariates.

