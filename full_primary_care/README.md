# full_primary_care

Contains scripts, QC log files, and results for analyses run on all phenotypes (commonly measured biomarkers, anthropometric, and spirometric measures) in primary care data; some combined with UK Biobank main phenotype files.  

## code_lists
There are 65 traits in primary care and 39 matched traits in UKB main phenotype and biomarker files. Code lists are curated from a combination of - manual, CALIBER, and Denaxas et al. 2020 lists - by Samvida Venkatesh and Kayesha Coley.

## scripts
1. **1_get_data.R** - Get primary care individual-level data for each of the 65 traits; perform initial QC to remove traits with fewer than 100 individuals with repeat measurements (hormones, anthropometric traits) or fewer than 1000 individuals with repeat measurements (all other traits). Returns 52 traits in a named list. 
2. **2_popn_QC.R** - Filter above data to only keep values that pass the age filter (ages 20-80yrs), filter out implausible values based on UKB main phenotypes stored in *code_lists/ukb_min_max.txt*, and extreme values +/- 5 S.D. away from mean value in each trait.

## qc
QC results from each step:

### popn_qc 
Separate files for each biomarker, detailing number of observations and number of individuals for each trait removed due to failing: age, plausibility, or extreme value filters (as detailed in *scripts/2_popn_QC.R*) 

