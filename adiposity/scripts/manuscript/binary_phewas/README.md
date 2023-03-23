Scripts to test the effect of rs429358 on 50+ longitudinal quantitative phenotypes.

1. **1_prep_longit_data.R** - Model linear slope change in quantitative traits with linear mixed effects models. Adjust the BLUP for confounders (baseline age, baseline age-squared, data provider, year of birth, and sex in sex-combined analyses). 
- **2_geno_sample_QC.R** As described in *../../GWAS/*.

Scripts in *longitudinal/*: Scripts to replicate obesity-change linear slope and cluster membership associations.
1. **1_prep_longit_data.R** - Get UK Biobank assessment centre values of BMI, weight, waist circumference (WC), and waist-hip-ratio (WHR), as well as self-reported weight change. Format for calculation of linear slope change and spline-based cluster membership (for BMI and weight).
2. **2_GWAS_sample_QC.R** - As described in *../../GWAS/*.
3. **3_genetic_assocn_linreg.R** - Test for the effect of rs429358 genotype, as calculated in *../post_GWAS/extract_dosages/*, on covariate-adjusted and appropriately transformed quantitative trait value slope BLUPs. 
4. **4_examine_plot_results.R** - Forest plots for effects.
