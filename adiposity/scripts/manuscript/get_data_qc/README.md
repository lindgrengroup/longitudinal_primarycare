Scripts in this folder:

1. **1_get_data.R** - Extract codes for all biomarkers of interest from full primary care data linked to UK Biobank, and perform initial QC to retain only biomarkers with > 1 measurement in > 1000 individuals in primary care (for longitudinal resource construction). 
2. **2_popn_QC.R** - Perform population-level QC to exclude records based on the following filters:
	- age not between 20-80 yrs
	- history or record of bariatric surgery (for obesity-related records only)
	- implausible values defined as outside +/-10% of the UKBB assessment centre minimum and maximum values
	- extreme values defined as > 5 SD away from population mean
3. **3_add_UKB_filter_longit.R** - Integrate UKBB assessment centre or biomarker measurements with linked primary care data (QC'd above) and filter to only retain individuals with longitudinal records (> 1 measurement of each biomarker).
4. **4_indiv-level_QC.R** - Remove records that are extreme on the individual level (i.e. records which cause large jumps in the time-series) by taking the log-fold-change between consecutive measurements, adjusted for time between consecutive measurements; any observations that are involved in an extreme jump > 3 SD away from population mean jump are removed.
5. **5_calculate_covariates.R** - Calculate trait-specific covariates, including: baseline age (age at first measurement), number of follow-up measures, length of follow-up, and baseline trait value. Store these as covariates for further analyses.

Code-lists for biomarkers are in the folder *code_lists/*, including a list of biomarkers in **qcd_traits_available.txt**, and their corresponding primary-care read codes (adapted from [Denaxas et al.](https://github.com/spiros/ukb-biomarker-phenotypes) and UKBB assessment centre codes