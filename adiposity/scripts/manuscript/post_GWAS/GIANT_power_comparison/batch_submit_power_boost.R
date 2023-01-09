# Author: Samvida S. Venkatesh
# Date: 15/04/2022

# Bulk submit power calculation jobs

compare_gp_giant <- read.table("compare_gp_giant.txt", sep = "\t", header = T, stringsAsFactors = F)

submission_script <- "scripts/submit_power_boost.sh"

for (i in 1:nrow(compare_gp_giant)) {
  job_options <- paste(
    "-v",
    paste0(
      "STRATA_NAME=\"", compare_gp_giant$strata[i], "\",",
      "META_FILENAME=\"", compare_gp_giant$meta_filename[i], "\",",
      "GP_N=\"", compare_gp_giant$gp_sample_size[i], "\",",
      "GP_FILENAME=\"", compare_gp_giant$gp_filename[i], "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}

