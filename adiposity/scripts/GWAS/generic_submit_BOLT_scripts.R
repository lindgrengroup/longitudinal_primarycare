# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit BOLT jobs
# Need to loop over all of the strata and parameter types in each strata

STRATA_NAMES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1
PARAMETERS <- c("lmm_intercepts", "lmm_slopes_adj_baseline", "cspline_intercepts")
# GWAS submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_helpers/1_perform_GWAS_BOLT.sh"
# Filtering submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_helpers/2_perform_BOLT_filtering.sh"

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/BOLT_helpers/2_perform_BOLT_filtering.sh"

for (s in STRATA_NAMES) {
  for (pr in PARAMETERS) {
    job_options <- paste(
      "-v",
      paste0(
        "STRATA=\"", s, "\",",
        "PARAMETER=\"", pr, "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
