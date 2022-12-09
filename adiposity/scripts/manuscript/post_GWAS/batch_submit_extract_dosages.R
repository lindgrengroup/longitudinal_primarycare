VARIDS <- c("rs429358", "rs8014554", "rs4420638", "rs707901")
CHRS <- c("19", "14", "19", "6")
submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/4_extract_dosages.sh"

for (i in 1:length(VARIDS)) {
  job_options <- paste(
    "-v",
    paste0(
      "VARID=\"", VARIDS[i], "\",",
      "CHR=\"", CHRS[i], "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
