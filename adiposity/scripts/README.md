Scripts in this folder:

## MAIN ANALYSES

### Cubic splines
1. **1_apply_time_splines.R** - Build cubic spline mixed effects models for effect of time on trait, allowing for fixed and random intercepts and spline effects of time. Run through models with increasing number of knots in the fixed effects (3-10) and random effects (1 - # fixed effects) to choose model with lowest BIC. Adjust for covariates (baseline age, baseline age-squared, and sex). Both sex-specific and sex-combined analyses. Calculate best linear unbiased predictor (BLUP) for each term as the fixed + random effect for each individual. 
2. **2_plot_cspline_predictions.R** - Plot model predictions for various subsets of individuals (random, high or low values of baseline trait, baseline age, covariates, etc.) 
3. **3_cBLUP_PCA_heatmaps.R** - Calculate PCs based on BLUPs from cubic spline models and plot heatmaps of BLUPs as well as PCA scree plots to visualise dimensionality reduction captured by PCs.
4. **4_*.R** - Try various clustering techniques: Gaussian mixture modelling, bootstrap number of clusters, and k-means clustering - all based on intercept and BLUPs from models.
5. **5_plot_clustering_results...** - Once clusters have been determined, plot modelled trajectories, observed trajectories, and associations between clusters and various covariates.
6. **6_remap_clusters.R** - If clusters have to be made to correspond across different strata.

### Cross-sectional (baseline)
1. **6_baseline_clusters.R** - Cluster the adjusted baseline trait value for each individual, stepping through 1:6 clusters allowing for varying shape, variance, etc. Assign individuals to clusters.


## PLOTTING AND DESCRIPTIVE STATISTICS

### Raw data - 
1. **report_characteristics.R** - Generate tables and figures for various descriptive characteristics of raw data.
2. **disease_associated_trajectories.R** - For a sample of diseases, plot raw data trajectories in individuals with and without disease.
### Modelling results - 
Change these scripts for each set of models, linear, cubic splines, etc.) - 
1. **plot_BLUP_distributions.R** - Plot randomly selected subset of BLUPs from models to check distributions as well as relationships between BLUPs and model covariates.
2. **plot_all_model_predictions.R** - Plot individual-level model predictions for randomly selected individuals, those with very few or many repeat measures, at the tails of model BLUPs, etc. for both linear and cubic spline models.
3. **trajectories_by_covariates/** - Plot adiposity trajectories (modelled by cubic splines) in each group of a covariate (ex. birth cohort, sex, disease chapter, etc.)
### Clustering results -  
1. **plot_clustering_results.R** - Once individuals have been assigned to a group/cluster belonging, plot various trajectories, refit models within each cluster, etc. 
2. **plot_cluster_disease_enrichment.R** - After testing for disease enrichment (see *../phenotype_enrichment/*), plot results in Miami-style plots.

## GWAS

Refer to README in */GWAS/*

## TRAJGWAS
Scripts to perform trajectory GWAS as outlined in https://github.com/OpenMendel/TrajGWAS.jl (Ko et al. 2022). Returns results for beta (mean), tau (variance), and joint effects of beta and tau for each SNP on trait. 
1. **1_sample_qc.R** - Prepare data and covariates in long format for TrajGWAS.
2. **2_fit_null_model.jl** - Use TrajGWAS Julia package to fit mixed model with fixed and random effects (including covariates and random intercept and random time-slope within each individual) but without any SNPs - null model. 
3. **3_spatest.jl** - Perform GWAS on each chromosome with the saddle point approximation implemented in Julia TrajGWAS pacakge.
4. **submit_4_plot_results.sh** and **4_plot_gwas_results.R** - Filter SNPs and plot QQ-plots and Manhattan plots.

## VARIANT SUBSET ASSOCIATIONS
Scripts to perform genetic associations with only a subset of variants known to associate with a range of obesity-traits (BMI, weight, waist circumference, waist-hip ratio, etc.)
1. **1_get_gwascat_obesity_associations.R** - Extract SNPs associated with obesity traits in GWAS catalog (downloaded 2021/11/02, build hg38). If missing chr and pos info, add this in. List of obesity traits defined in associated file **gwascat_obesity_traits.txt**. Write associations to bed file to lift over to hg19.
2. **2_liftover_hg38_to_hg19.sh** - Convert hg38 chr/pos information SNPs to hg19.
3. **3_collate_snp_list_per_chr.R** - Write list of rsids and variant ids per chromosome to extract from UKB bgen files.

## Post-GWAS
1. **1_perform_finemapping.sh** - Uses Duncan's pipeline (see here: https://github.com/astheeggeggs/pipeline) to finemap causal variants in loci from filtered GWAS summary statistics. 
2. **2_plot_locuszoom.sh** - Following finemapping, provide LD matrix, summary stats, and finemapped variants to local version of locuszoom for loci plots. This can be extended to highlight independently associated variants. 
3. **3_calculate_heritability.sh** - Get overall SNP-based heritability of phenotype using LDSC from summary statistics. 
4. **4_extract_dosages.sh and 4_submit_extract_dosages.R** - Get dosages (0,1,2) to find homozygous REF/ALT and heterozygous individuals at lead SNPs of interest. 

Other helpers:
**8_calculate_power_boost** - Compare chi-square statistics from linear mixed model intercept GWASs to GIANT GWASs to determine increase in effective sample size.
**classify/characterise_reported_novel_variants.R** - Classify variants as "reported" (previously published), "refined" (correlated with previously published variants, larger effect size than previously published correlated variants, and retains effect when conditioned on previously published correlated variants), or "novel" (uncorrelated with previously published variants, and conditionally independent of any previously published variants in the region). Characterise with global MAF, nearest gene, and variant consequence.

## PGS ASSOCIATIONS (in collaboration with Frederik Heymann Lassen)
See here for generation of PGS - https://github.com/frhl/wes_ko_ukbb/tree/speedy_speedos/scripts/prs
1. **plot_PGS_in_clusters.R** - Compare polygenic score distributions for various traits across the clusters.

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

