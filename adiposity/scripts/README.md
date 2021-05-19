Scripts in this folder:

## MAIN ANALYSES

### General
1. **1_get_phenotypes.R** - Gather primary care trait measurements and add same traits from UKBIOBANK, removes individuals without repeat measures.
2. **2_clean_phenotypes.R** - For each trait, QC on age, implausible and extreme values, and remove individuals without repeat measures. Flags pregnancy. Inter-convert all BMI and weight measures.
3. **3_plot_QC.R** - Visual inspection for unrealistically large jumps between time-points. Remove measurements causing unrealistic jump.
4. **4_calculate_covariates.R** - Calculate individual-level phenotyping (baseline age, baseline BMI, follow-up years, etc.) and genotyping (array, PCs) covariates. Stratify on sex and ancestry.

### Linear mixed models
5. **5_raw_slopes.R** - Perform mixed effects regression with fixed and random effects for individual and age, regressing adiposity trait on age. Calculate BLUP for raw slope for each individual, remove outliers > 5 S.D. away from mean in each stratum. 
6. **6_adjust_slopes.R** - Run linear models to adjust raw BLUPs for a range of covariates with nested models. Calculate variance explained by baseline model covariates and compare models with ANOVA. (6a) picks the best model to carry forward.
7. **7_adjust_slope_groups.R** - Calculate groups (gainers, quartiles, etc.).

### Spline mixed models
5. **5_raw_slopes.R** - Perform mixed effects regression with splines (NKNOTS = 8 for BMI and 3 for WHR) with fixed and random effects for individual and age, regressing adiposity trait on age. Calculate BLUP for raw slope for each individual, flag outliers > 5 S.D. away from mean in each stratum. (5a) plots the predictions from these models.
6. **6_adjust_slopes.R** - Run linear models to adjust raw BLUPs for covariates picked above (baseline age, baseline age-squared, baseline adiposity trait, number of follow-up years and number of follow-up measures, first 21 genetic principal components, and sex in sex-combined analyses).

## PLOTTING, DESCRIPTIVE TABLES, ETC.
1. **compare_adj_slope_models.R** - ANOVA, AIC, etc. to compare models for BLUP adjustment, plot mean trajectories by different model quartiles to compare effect of adjustment.
2. **report_characteristics.R** - Generate tables and figures for various descriptive characteristics of raw data, raw BLUPs (and trajectories), and adjusted BLUPs (and trajectories).

## PHENOTYPE ENRICHMENT
1. **female_phenotype_enrichment.R** - Get age at menarche, age at menopause, and number of live births from UKB to compare distributions in different BLUP quartiles in women.
2. **phenotype_enrichment.R** - Fisher's enrichment test for disease codes (Spiros, primary + secondary care) in specific groups (ex. gainers) as compared to cohort.

## GWAS
1. **prepare_GWAS_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Get GWAS phenotype - RINTed adjusted BLUP - for QCd individuals in full cohort as well as in only gainers.
2. **perform_GWAS_{BOLT/SUGEN}.sh** - Perform genome-wide association study using specified software.

## MULTIVARIATE ANALYSES, PLOTS, ENRICHMENT, ETC.
1. **focus_BMI_WHR_overlap.R** - Data summaries and enrichment tests for individuals with both BMI and WHR trajectories in the data.
2. **multivariate_report_characteristics.R** - Data summaries, overlaps, and trajectories for individuals who have more than one adiposity trait measured in the data.

## ARCHIVED:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Removes population-level noise by calculating thresholds based on UKBIOBANK-measured BMI and trajectories of extreme individuals 
3. **bmi_post_stage1_stats.R** - Plots basic information such as age distribution, sex distribution, longitudinal traits, etc. of cleaned BMI measures; same script used after individual outlier removal
4. **bmi_individual_outlier_removal.R** - Removes individual-level noise by calculating sequential logFC (fold-change between consecutive measurements) as well as accounting for time between consecutive measurements; individuals in the tails of these distributions are highlighted and the noisy observation (farthest from individual median or UKBIOBANK BMI) is removed
5. **bmi_trajectories.R** - Plots trajectories for individuals stratified by number of observations, overall change in BMI, etc. 
