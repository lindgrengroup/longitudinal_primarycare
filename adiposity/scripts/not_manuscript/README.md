## The common genetic architecture of adult obesity change

Code for analyses in article here - 

### Cite as: 

## Structure
Scripts are organised by directories corresponding to sections in the manuscript.

- *get_data_qc/* - Identification and quality control of longitudinal records
- *linear_mixed_models/* - Linear mixed-effects models to define baseline obesity (intercept) and obesity-change (slope)
- *high_dim_splines/* - Regularised spline models for non-linear obesity trajectories, and soft clustering of individuals
- *GWAS/* - Genome-wide association studies and finemapping
- *post_GWAS/* - Analyses that require summary statistics from GWAS, such as power comparison to previous studies, heritability and genetic correlation calculations, sex-heterogeneity testing, etc. 
- *replication_ukb/* - Replication of models and genome-wide significant associations in UK Biobank hold-out sets
- *longit_phewas/* - Longitudinal phenome-wide association studies for rs429358
- *non_wb_ancestry/* - Identification and quality control of longitudinal records from individuals not of white British ancestry, replication of models and genome-wide significant associations in this subset
- *manuscript_figures/* - Scripts to generate figures and supplementary figures in manuscript
- *manuscript_tables/* - Scripts to generate tables in manuscript, if not already in other directories


## MAIN ANALYSES

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

