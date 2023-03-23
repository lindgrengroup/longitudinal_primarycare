Scripts in this folder:

Following reviewer request, comparison of LMM intercepts to average trait across all measurements.

1. **get_average_compare_lmm_intercept.R** - Calculate average from across measurements and plot phenotype average vs LMM BLUP intercept.
2. **2_RINT_average_trait.R** - Adjust for covariates and rank-based inverse normal transform the residuals to carry forward for GWAS.

3. **3_extract_gws_snps.sh** - Post-GWAS and filtering, curate a list of SNPs that are genome-wide significant (P < 5E-8) in either the LMM intercept GWAS or the average trait GWAS - extract summary statistics for these SNPs.
4. **4_plot_scatter_compare_effects.R** - Plot scatter of betas from the two GWASs to compare effect sizes between LMM intercept GWAS and average trait GWAS.

