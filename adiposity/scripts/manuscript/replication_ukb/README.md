Scripts to replicate the genome-wide significant associations discovered in the primary-care + UKBB data in UKBB assessment centre data alone.

1. **debias_beta_check_replication_sample_sizes.R** - Estimate the required sample sizes to replicate linear slope associations after adjusting the discovery effect size for winner's curse. Debiasing script adapted from [Palmer *et al.* (2017)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006916).

Scripts in *cross_sectional/*: Scripts to replicate baseline obesity associations.

1. **1_prep_cross_sec.R** - Get UK Biobank cross-sectional value by taking the assessment centre visit measurement of BMI or weight, or if there are multiple measurements, the median measurement. Adjust for confounders (age, age-squared, data provider, year of birth, and sex in sex-combined analyses), and rank-based inverse normally transform for GWAS. 
2. **2_sample_QC.R** and **3_perform_GWAS_BOLT.sh** - Perform sample-level quality control for the UKB held-out samples for GWAS. Perform GWAS with BOLT-LMM. Similar parameters as detailed in *../../GWAS/*

Scripts in *longitudinal/*: Scripts to replicate obesity-change linear slope and cluster membership associations.
1. **1_prep_longit_data.R** - Get UK Biobank assessment centre values of BMI, weight, waist circumference (WC), and waist-hip-ratio (WHR), as well as self-reported weight change. Format for calculation of linear slope change and spline-based cluster membership (for BMI and weight).
2. **2_GWAS_sample_QC.R** - As described in *../../GWAS/*.

Calculate adiposity-change phenotypes (linear slope and posterior probability of cluster membership) following the scripts in */lmm_slopes/* and *highdim_splines/*. All these steps are outlined in more detail in the README files in *../../linear_mixed_models/* and *../../high_dim_splines/* where the discovery analysis scripts are stored.

**subset_replication_genotypes.sh** - Extract the list of SNPs associated with obesity-change trait in discovery analyses and genotypes/dosages at these SNPs for all individuals. 

3. **3_GWAS_trait_curation.R** - Gather the adiposity-change traits calculated above and adjust for confounders and RINT or transform posterior probability of cluster memberships. 
4. **4_variant_association_analysis.sh** - Since we are only testing the association of a subset of individuals with a subset of variants, use PLINK to perform association analysis, adjusting for appropriate covariates. Also run **longit_effect_main_dat.R**, which tests for the effect of genotypes of interest on waist-circumference-change, waist-to-hip-ratio-change, and self-reported weight change, and **4_rs429358_no_dementia.R**, which tests for the association of this single SNP with all adiposity change phenotypes, including self-reported weight change, after excluding individuals with dementia (as calculated in *../../get_data_qc/*)
