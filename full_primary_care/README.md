# full_primary_care

Contains scripts, QC log files, and results for analyses run on all phenotypes (commonly measured biomarkers, anthropometric, and spirometric measures) in primary care data; some combined with UK Biobank main phenotype files.  

## code_lists
There are 65 traits in primary care and 37 matched traits in UKB main phenotype and biomarker files. Code lists are curated from a combination of - manual, CALIBER, and Denaxas et al. 2020 lists - by Samvida Venkatesh and Kayesha Coley.

## scripts
1. **1_get_data.R** - Get primary care individual-level data for each of the 65 traits; perform initial QC to remove traits with fewer than 100 individuals with repeat measurements (hormones, anthropometric traits) or fewer than 1000 individuals with repeat measurements (all other traits). Returns 52 traits in a named list. 
2. **2_popn_QC.R** - Filter above data to only keep values that pass the age filter (ages 20-80yrs), filter out implausible values based on UKB main phenotypes stored in *code_lists/ukb_min_max.txt*, and extreme values +/- 5 S.D. away from mean value in each trait.
3. **3a_ and 3b_** - Calculate descriptive factors on baseline primary care data, such as age at first measurement, median baseline trait value, etc. Perform test for heterogeneity between sexes. 
4. **4_add_UKB_filter_longit.R** - Add UK Biobank measurements from assessment centres (where available), to individuals with primary care recordings and long-format data. Remove any individuals that have < 2 measurements. 
5. **5_indiv_level_QC.R** - Longitudinal filtering on individual measurements, removing timepoints where there are > 1 measurements at the same timepoint; also check for extreme "jumps", i.e. log fold-change between consecutive measurements in the data, to exclude timepoints with extreme changes.
6. **6_calculate_covariates.R** - Get trait-specific covariates, such as baseline age (age at first measurement of trait), baseline trait value, number of follow-up measurements, number of years between first and last measurement, etc.
7. **plot_age_at_diagnosis.R** - Across all phenotypes, plot distribution of age at diagnosis for males and females.
7. **plot_trajectories.R** - For a list of IDs with classification and biomarkers passed as arguments, this script plots (a) distribution of covariates such as baseline age, baseline biomarker, # follow-ups, etc. in each class, (b) randomly sampled individual trajectories of specified biomarkers in each class, and (c) group-level mean trajectories of specified biomarkers in each class.  

## qc
QC results from each step:

### popn_qc 
Separate files for each biomarker, detailing number of observations and number of individuals for each trait removed due to failing: age, plausibility, or extreme value filters (as detailed in *scripts/2_popn_QC.R*) 

