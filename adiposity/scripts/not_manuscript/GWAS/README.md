Scripts in this folder:

1. **1_sample_QC.R** - Genotyping-related sample quality control (ex. remove individuals with reported and genotyped sex mismatches, retain only individuals in the white British ancestry subset, samples with poor heterozygosity or missingness, etc.). Save sample ids that pass QC along with genotyping-related covariates, i.e. genotyping array and UKB assessment centre. This only needs to be run once for each strata, as all models share the same individuals.
2. Scripts to create GWAS outcomes -
**2_create_hard_cluster_traits.R** - Assign individuals to cluster with greatest posterior probability after soft clustering --> feed this foward to ordinal logistic regression GWAS with POLMM. \
**2_create_soft_cluster_traits.R** - Calculate logit-transformed squeezed probabilities of belonging to k1, k1 or k2, etc. to create traits suitable for linear regression GWAS with BOLT-LMM. \
**2_RINT_lmm_traits.R** - Adjust the linear mixed model intercept and slope BLUPs for covariates and rank-based inverse normal transform the residuals within each sex strata to create traits for linear regression GWAS with BOLT-LMM. \
3. **generic_submit_BOLT_scripts.R** - Loops over all strata and parameters (lmm intercepts, lmm slopes adjusted for baseline, cubic spline intercepts) to submit BOLT GWAS, filtering, etc. Relies on helper scripts in *./BOLT_helpers/*. 
4. **softprob_linear_regression_GWAS.sh** - Follows the same structure as *./BOLT_helpers/1_perform_GWAS_BOLT.sh* but includes the covariates needed for soft-probability GWAS as these traits are not adjusted for covariates prior to GWAS. 

Scripts in *./BOLT_helpers/*:

1. **1_perform_GWAS_BOLT.sh** - Take in arguments passed by BOLT submission script (for strata and parameter) to perform GWAS under the linear mixed model framework in BOLT, using imputed genotypes from UK Biobank. 
2. **2_perform_BOLT_filtering.sh** - Shell script to (1) Perform initial filtering by MAF > 1%, HWE p-value > 1E-6, INFO > 0.8, missingness < 5%, and bi-allelic SNPs, and (2) execute the filtering wrapper script **2_BOLT_filtering_wrapper.R**. Output: filtered gzipped GWAS results in a format suitable for FUMA, and plot QQ plots and Manhattan plots.
3. **3_perform_finemapping_BOLT.sh** - Calls Duncan's finemapping pipeline (using the FINEMAP software) to identify putative causal SNPs from GWAS summary statistics. 

Scripts in *./POLMM_helpers/*:

1. **1_POLMM_step1_wrapper.R** - Called by the shell script **1_perform_POLMM_step1.sh** to prepare for single-variant testing in an ordinal logistic regression framework by specifying sparse GRM, genotype files, etc.
1. **2_POLMM_step2_wrapper.R** - Called by the shell script **2_perform_POLMM_step2.sh** to perform single-variant testing in an ordinal logistic regression framework.
2. **3_perform_POLMM_filtering.sh** - Shell script to (1) Perform initial filtering by MAF > 1%, HWE p-value > 1E-6, INFO > 0.8, missingness < 5%, and bi-allelic SNPs, and (2) execute the filtering wrapper script **3_POLMM_filtering_wrapper.R**. Output: filtered gzipped GWAS results in a format suitable for FUMA, and plot QQ plots and Manhattan plots.
3. **4_perform_finemapping_POLMM.sh** - Calls Duncan's finemapping pipeline (using the FINEMAP software) to identify putative causal SNPs from GWAS summary statistics. 
