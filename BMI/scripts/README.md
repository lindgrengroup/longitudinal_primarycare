Scripts in this folder to be run in the following order:

1. **bmi_basic_stats.R** - Gathers primary care BMI measures, adds information from UKBIOBANK on sex, date of birth, and mean UKBIOBANK BMI, calculates age at measurement; then plots basic information such as age distribution, sex distribution, etc. of primary care BMI measures
2. **bmi_popn_outlier_removal.R** - Remove errors in primary care BMI measures based on upper and lower thresholds calculated as +/- 10% of the most extreme values seen in UKBIOBANK BMI data; modify these thresholds based on trajectories of the most extreme 0.1% of individuals in UKBIOBANK BMI; plot error characteristics
