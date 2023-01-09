# Author: Samvida S. Venkatesh
# Date: 01/03/2022

# R wrapper to submit qctool jobs
# Loop over all lead variants in each parameter (across all strata)

submission_script <- "scripts/extract_dosages.sh"

lead_variants <- read.table("lead_snps_to_replicate.txt",
                            sep = "\t", header = T, stringsAsFactors = F)

for (vi in 1:nrow(lead_variants)) {
  snp <- lead_variants$SNP[vi]
  if (grepl("^chr", snp)) {
    snp <- gsub("chr", "", snp)
    snp <- paste0(snp, "_", lead_variants$Tested_Allele[vi], "_", lead_variants$Other_Allele[vi])
  }
  
  job_options <- paste(
    "-v",
    paste0(
      "VARID=\"", snp, "\",",
      "CHR=\"", lead_variants$CHR[vi], "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
