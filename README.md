# longitudinal_primarycare

## code_lists
Contains primary care (read V2, read V3) and UKB main phenotype file fields for commonly measured biomarkers, anthropometric measures, and spirometric measures. There are 65 traits in primary care and 21 (18 biomarkers + 3 anthropometric measures) matched traits in UKB main phenotype files. Code lists are curated from a combination of - manual, CALIBER, and Denaxas et al. 2020 lists - by Samvida Venkatesh and Kayesha Coley. 

## general_scripts
General scripts:
1. **bariatric_surgery_annotation.R** - Produces list of individuals with bariatric surgery records in UKB or primary care and date of surgery if known.
2. **calculate_demographic_covariates.R** - Annotates individual IDs with sex (confirmed by UKB QC file), ancestry, median height, and first 21 genetic PCs from UKB. 
3. **gp_age_sex_BMI_annotation.R** - Annotates the entire GP clinical UKBIOBANK file with UKBIOBANK-phenotype details on sex, mean BMI from UKBIOBANK measurements, and DOB (with age calculation).
4. **pregnancy_annotation.R** - Produces list of individuals with pregnancy records in UKB or primary care and estimated dates of pregnancy.
5. **spiros_phenotype_codes_cleaning.R** - Adds ICD9 codes to the ICD10 codes and priamry care read CTV3 codes to the V2 codes in the cleaned phenotype list produced by Spiros Denaxas group (https://github.com/spiros/chronological-map-phenotypes).
6. **build_eid_phenotype_matrix.R** - Constructs matrix of 0/1 for presence of disease phenotype code (cleaned, Spiros) in each individual in UKB secondary data and primary care data.
7. **time_to_event_from_**.R** - Calculate age at first diagnosis for each disease phenotype code (cleaned, Spiros) in each individual in UKB primary care data and HES data.
8. **ukb_extract_wb_ids.R** - Get subset of individuals of White British ancestry (as genetically identified by Bycroft et al. 2018).
9. **hes_get_age_event.R** - Annotate HES records in UK Biobank with age at event, calculated with UKBIOBANK-based DOB, choosing the earliest of multiple dates for the same record.

## scripts
Remaining scripts organised by trait:
## multivariate
## adiposity
## hormones
