# longitudinal_primarycare
Scripts for assessing longitudinal quantitative traits in UKBIOBANK-linked primary care data

General scripts in /general_scripts:
1. **bariatric_surgery_annotation.R** - Produces list of individuals with bariatric surgery records in UKB or primary care and date of surgery if known.
2. **calculate_demographic_covariates.R** - Annotates individual IDs with sex (confirmed by UKB QC file), ancestry, median height, and first 21 genetic PCs from UKB. 
3. **gp_age_sex_BMI_annotation.R** - Annotates the entire GP clinical UKBIOBANK file with UKBIOBANK-phenotype details on sex, mean BMI from UKBIOBANK measurements, and DOB (with age calculation).
4. **pregnancy_annotation.R** - Produces list of individuals with pregnancy records in UKB or primary care and estimated dates of pregnancy.

Remaining scripts organised by trait /scripts:

## multivariate
## adiposity
## hormones
