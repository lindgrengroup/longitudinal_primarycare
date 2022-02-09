# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit SAIGE jobs
# Need to loop over all of the strata and clusters in each strata

STRATA_NAMES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1
CLUSTERS <- 1:6
submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/perform_SAIGE_step1.sh"

for (s in STRATA_NAMES) {
  for (c in CLUSTERS) {
    job_options <- paste(
      "-v",
      paste0(
        "STRATA=\"", s, "\",",
        "KI=\"", c, "\""
      )
    )
    job_submission <- paste("qsub", job_options, submission_script)
    system(job_submission)
    print(job_submission)
  }
}
