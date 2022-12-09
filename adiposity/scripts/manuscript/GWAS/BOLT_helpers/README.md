Scripts in this folder:

**1_perform_GWAS_BOLT.sh** - Take in arguments passed by BOLT submission script (for strata and parameter) to perform GWAS under the linear mixed model framework in BOLT, using imputed genotypes from UK Biobank. Perform initial filtering by MAF > 1%, HWE p-value > 1E-6, INFO > 0.8, missingness < 5%, and bi-allelic SNPs.

**2_perform_BOLT_filtering.sh** - Shell script to (1) filter and log filtering GWAS results from BOLT, and (2) execute the filtering wrapper script below with command-line arguments submitted by the submit scripts in *../* that loop over different strata and parameters. Output: filtered gzipped GWAS results. \
**2_BOLT_filtering_wrapper.R** - Script to perform QC on *filtered gzipped* GWAS results, output summary statistics in a format suitable for FUMA, and plot QQ plots and Manhattan plots. Parses arguments provided from the command line in previous script.
