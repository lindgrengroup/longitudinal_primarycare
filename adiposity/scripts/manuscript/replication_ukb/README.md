Scripts to replicate the genome-wide significant associations discovered in the primary-care + UKBB data in UKBB assessment centre data alone.

Scripts in *cross_sectional/*: Scripts to replicate baseline obesity associations.

- **1_prep_cross_sec.R** - Get UK Biobank cross-sectional value by taking the assessment centre visit measurement of BMI or weight, or if there are multiple measurements, the median measurement. Adjust for confounders (age, age-squared, data provider, year of birth, and sex in sex-combined analyses), and rank-based inverse normally transform for GWAS. 
- **2_sample_QC.R** and **3_perform_GWAS_BOLT.sh** - Perform sample-level quality control for the UKB held-out samples for GWAS. Perform GWAS with BOLT-LMM. Similar parameters as detailed in *../../GWAS/*




- **gcta_cojo_conditional_analysis.sh** - For a given list of SNPs (containing a SNP of interest and all obesity-associated SNPs within 500kb), create LD matrix based on individuals of White British ancestry to perform conditional analysis with GCTA-COJO. All collinear SNPs are added to the "reported" SNP list; if the SNP of interest remains as an independent association, it may be refined or novel -- calculate conditional effect of all independent SNPs in the region.
- **characterise_refined_novel_variants.R** - From the conditional effects calculated above, classify SNPs as reported, refined, or novel based on the "significantly stronger effect" and "conditionally independent" criteria using t-tests. Access "buddy" SNPs for each refined and novel SNP, i.e. most highly correlated or nearest obesity-associated variant.

Scripts in *GIANT_power_comparison/*: Scripts to assess ratio of chi-squared statistics from in-house BMI intercept GWAS to published [GIANT 2019 meta-analysis of BMI](https://academic.oup.com/hmg/article/28/1/166/5098227), to calculate effective sample size boost from repeat measurements.  Calculation performed in **calculate_power_boost_chisq.R**, summary statistics formatted and job submitted in **submit_power_boost.sh**, via batch submission script **batch_submit_power_boost.R**.

Scripts in *LDSC_r2_hg/*: Scripts to calculate SNP-based heritability and genetic correlation between: obesity-intercept and obesity-change traits, as well as between BMI and weight for the different change phenotypes. Munge sumstats into the correct format, and use SNPs in the pan-UKBB LD panel (EUR) for 1 million HapMap variants.

Scripts in *extract_dosages/*: Scripts to extract genotype dosages at variants of interest (lead SNPs) using QCTOOLS.
