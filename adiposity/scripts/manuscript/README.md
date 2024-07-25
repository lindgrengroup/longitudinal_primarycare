# The genetic architecture of changes in adiposity during adulthood

Code for analyses in published article here - https://www.nature.com/articles/s41467-024-49998-0

Cite as: Venkatesh, S.S., Ganjgahi, H., Palmer, D.S. et al. Characterising the genetic architecture of changes in adiposity during adulthood using electronic health records. Nat Commun 15, 5801 (2024). https://doi.org/10.1038/s41467-024-49998-0

## Structure
Scripts are organised by directories corresponding to sections in the manuscript.

- *get_data_qc/* - Identification and quality control of longitudinal records
- *linear_mixed_models/* - Linear mixed-effects models to define baseline obesity (intercept) and obesity-change (slope)
- *high_dim_splines/* - Regularised spline models for non-linear obesity trajectories, and soft clustering of individuals
- *GWAS/* - Genome-wide association studies and finemapping
- *post_GWAS/* - Analyses that require summary statistics from GWAS, such as classification of variants as reported or novel, power comparison to previous studies, heritability and genetic correlation calculations, and extracting genotypes/dosages of lead SNPs. 
- *replication_ukb/* - Replication of models and genome-wide significant associations in UK Biobank hold-out sets
- *longit_phewas/* - Longitudinal phenome-wide association studies for all lead variants
- *non_wb_ancestry/* - Identification and quality control of longitudinal records from individuals not of white British ancestry, replication of models and genome-wide significant associations in this subset
- *manuscript_figures/* - Scripts to generate figures and supplementary figures in manuscript, if not already in other directories

## To address review, the following directories were added
- *TrajGWAS/* - Using the [TrajGWAS package](https://github.com/OpenMendel/TrajGWAS.jl) for GWAS of intra-individual mean and variance
- *binary_phewas/* - Binary disease phenome-wide association studies for all lead variants
- *sensitivity_compare_FU_adjustment/* - To assess the difference between follow-up adjusted and unadjusted GWASs
- *sensitivity_compare_intercept_average/* - To assess the difference between the LMM intercept and average BMI GWASs
- *sensitivity_longevity_snps/* - Testing whether lifespan-associated variants are also associated with adiposity-change
- *sensitivity_rs429358/* - Additional checks for survivorship bias for this lifespan-associated variant
