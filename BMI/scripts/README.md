Scripts in this folder:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Removes population-level noise by calculating thresholds based on UKBIOBANK-measured BMI and trajectories of extreme individuals 
3. **bmi_post_stage1_stats.R** - Plots basic information such as age distribution, sex distribution, longitudinal traits, etc. of cleaned BMI measures; same script used after individual outlier removal
4. **bmi_individual_outlier_removal.R** - Removes individual-level noise by calculating sequential logFC (fold-change between consecutive measurements) as well as accounting for time between consecutive measurements; individuals in the tails of these distributions are highlighted and the noisy observation (farthest from individual median or UKBIOBANK BMI) is removed
5. **bmi_trajectories.R** - Plots trajectories for individuals stratified by number of observations, overall change in BMI, etc. ADD MORE AS MORE WAYS TO LOOK AT DATA
