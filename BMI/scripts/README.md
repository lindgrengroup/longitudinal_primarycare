Scripts in this folder to be run in the following order:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Removes population-level noise by calculating thresholds based on UKBIOBANK-measured BMI and trajectories of extreme individuals 
3. **bmi_post_stage1_stats.R** - Plots basic information such as age distribution, sex distribution, longitudinal traits, etc. of cleaned BMI measures
4. **bmi_individual_outlier_removal.R** - IN PROGRESS; Method 1 - min-max logFC
