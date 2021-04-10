Scripts in this folder:

1. **1_get_phenotypes.R** - Gather primary care trait measurements and adds same traits from UKBIOBANK, removes individuals without repeat measures.
2. **2_clean_phenotypes.R** - For each trait, QC on age, implausible and extreme values, and remove individuals without repeat measures. Flags pregnancy. Inter-converts all BMI and weight measures.
3. **3_plot_QC.R** - Visual inspection for unrealistically large jumps between time-points. Remove measurements causing unrealistic jump.
4. **4_calculate_covariates.R** - Calculate individual-level phenotyping (baseline age, baseline BMI, follow-up years, etc.) and genotyping (array, PCs) covariates. Stratify on sex and ancestry.
5. **5_raw_slopes.R** - Perform mixed effects regression with fixed and random effects for individual and age, regressing adiposity trait on age. Calculate BLUP for raw slope for each individual, remove outliers > 5 S.D. away from mean in each stratum. 
6. **6_adjust_slopes.R** - Run linear models to adjust raw slopes for a range of covariates with nested models. Calculate variance explained by baseline model covariates and compare models with ANOVA.
7. **7_adjust_slope_groups.R** - Choose best model and calculate groups (gainers, quartiles, etc.).
8. **phenotype_enrichment.R** - Fisher's enrichment test for disease codes (Spiros, primary + secondary care) in specific groups (ex. gainers) as compared to cohort.
9. **prepare_GWAS.R** - Perform genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Get GWAS phenotype - RINTed adjusted slope - for QCd individuals in full cohort as well as in only gainers.
10. **report_characteristics.R** - Generate tables and figures for various descriptive characteristics of raw data, raw slopes (and trajectories), and adjusted slopes (and trajectories).
11. **focus_BMI_WHR_overlap.R** - Data summaries and enrichment tests for individuals with both BMI and WHR trajectories in the data.

ARCHIVED:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Removes population-level noise by calculating thresholds based on UKBIOBANK-measured BMI and trajectories of extreme individuals 
3. **bmi_post_stage1_stats.R** - Plots basic information such as age distribution, sex distribution, longitudinal traits, etc. of cleaned BMI measures; same script used after individual outlier removal
4. **bmi_individual_outlier_removal.R** - Removes individual-level noise by calculating sequential logFC (fold-change between consecutive measurements) as well as accounting for time between consecutive measurements; individuals in the tails of these distributions are highlighted and the noisy observation (farthest from individual median or UKBIOBANK BMI) is removed
5. **bmi_trajectories.R** - Plots trajectories for individuals stratified by number of observations, overall change in BMI, etc. ADD MORE AS MORE WAYS TO LOOK AT DATA
