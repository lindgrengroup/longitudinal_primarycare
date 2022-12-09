Scripts in this folder:

1. **1_get_data.R** - Extract codes for all biomarkers of interest from full primary care data linked to UK Biobank, and perform initial QC to retain only biomarkers with > 1 measurement in > 1000 individuals in primary care (for longitudinal resource construction). 
2. **2_popn_QC.R** - Perform population-level QC to exclude records based on the following filters:
	- age not between 20-80 yrs
	- history or record of bariatric surgery (for obesity-related records only)
	- implausible values defined as outside +/-10% of the UKBB assessment centre minimum and maximum values
	- extreme values defined as > 5 SD away from population mean

Code-lists for biomarkers are in the folder *code_lists/*, including a list of biomarkers in **qcd_traits_available.txt**, and their corresponding primary-care read codes (adapted from [Denaxas et al.](https://github.com/spiros/ukb-biomarker-phenotypes)) and UKBB assessment centre codes
