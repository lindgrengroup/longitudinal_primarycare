# Author: Samvida S. Venkatesh
# Date: 08/02/2022

# R wrapper to submit clustering jobs
# Need to loop over all of the strata and parameter types in each strata

strata_db <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/clustering/mclust_parameters.txt",
                        sep = "\t", header = T, stringsAsFactors = F)

script_to_sub <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/submit_gmm_soft_clustering.sh"


for (i in 1:nrow(strata_db)) {
  job_options <- paste(
    "-v",
    paste0(
      "STRATA=\"", strata_db$strata[i], "\",",
      "NCLUST=\"", strata_db$nclust[i], "\""
    )
  )
  job_submission <- paste("qsub", job_options, script_to_sub)
  system(job_submission)
  print(job_submission)
}

