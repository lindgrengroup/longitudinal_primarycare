Scripts in this folder:

**1/2_SAIGE_step1/2_wrapper.R** - Script to call models for SAIGE step 1 and step 2, with the functions fitNULLGLMM and SPAGMMATtest respectively, by parsing arguments provided from the command line and concatenating file paths. \
**1/2_perform_SAIGE_step1/2.sh** - Shell script to execute the wrapper scripts above with command-line arguments submitted by the submit scripts in *../* that loop over different strata and clusters.

**3_perform_SAIGE_filtering.sh** - Shell script to (1) concatenate, filter, and log filtering GWAS results from SAIGE step 2, and (2) execute the filtering wrapper script below with command-line arguments submitted by the submit scripts in *../* that loop over different strata and clusters. Output: filtered gzipped GWAS results. \
**3_SAIGE_filtering_wrapper.R** - Script to perform QC on *filtered gzipped* GWAS results, output summary statistics in a format suitable for FUMA, and plot QQ plots and Manhattan plots. Parses arguments provided from the command line in previous script.
