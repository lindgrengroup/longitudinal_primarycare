# Author: Samvida S. Venkatesh
# Date: 08/02/2022
# Adapted from: Duncan Palmer (https://github.com/astheeggeggs/SAIGE_gene_munging)

# R wrapper to submit SAIGE jobs
# Need to loop over all of the strata and clusters in each strata

STRATA_NAMES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/strata_list.txt")$V1

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/scripts/submit_3_spatest.sh"

for (s in STRATA_NAMES) {
  job_options <- paste(
    "-v",
    paste0(
      "STRATA=\"", s, "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
