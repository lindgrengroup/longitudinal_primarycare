VARIDS <- "rs429358"
CHRS <- "19"

submission_script <- "scripts/extract_dosages.sh"

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
