# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit SAIGE jobs
# Need to loop over all of the strata and clusters in each strata

STRATA_NAMES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1
CLUSTERS <- 1:6
# Step 1 submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/1_perform_SAIGE_step1.sh"
# Step 2 submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/2_perform_SAIGE_step2.sh"
# Filtering submission script
# submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/3_perform_SAIGE_filtering.sh"

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/SAIGE_helpers/3_perform_SAIGE_filtering.sh"

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
