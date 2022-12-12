# The common genetic architecture of adult obesity change

Code for analyses in article here - 

Cite as: 

## Structure
Scripts are organised by directories corresponding to sections in the manuscript.

- *get_data_qc/* - Identification and quality control of longitudinal records
- *linear_mixed_models/* - Linear mixed-effects models to define baseline obesity (intercept) and obesity-change (slope)
- *high_dim_splines/* - Regularised spline models for non-linear obesity trajectories, and soft clustering of individuals
- *GWAS/* - Genome-wide association studies and finemapping
- *post_GWAS/* - Analyses that require summary statistics from GWAS, such as power comparison to previous studies, heritability and genetic correlation calculations, sex-heterogeneity testing, etc. 
- *replication_ukb/* - Replication of models and genome-wide significant associations in UK Biobank hold-out sets
- *longit_phewas/* - Longitudinal phenome-wide association studies for rs429358
- *non_wb_ancestry/* - Identification and quality control of longitudinal records from individuals not of white British ancestry, replication of models and genome-wide significant associations in this subset
- *manuscript_figures/* - Scripts to generate figures and supplementary figures in manuscript
- *manuscript_tables/* - Scripts to generate tables and supplementary tables in manuscript, if not already in other directories
