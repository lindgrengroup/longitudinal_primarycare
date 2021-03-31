# Author: Samvida S. Venkatesh
# Date: 30/03/21

library(tidyverse)

# Read data ----

# Read dictionary to get to phenotype lists
dictionary <- read.csv("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/dictionary.csv",
                         header = T)
PHENOTYPES <- dictionary$phenotype

# Read merged V2-V3 file
merged_v2v3 <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/merged_v2v3_codes.txt",
                          sep = "\t", header = T, 
                          quote = "", fill = F, comment.char = "~")

# Read merged ICD9-ICD10 file
merged_icd <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/merged_icd9_10_codes.txt",
                         sep = "\t", header = T, 
                         quote = "", fill = F, comment.char = "~")

# Create new column for dictionary with annotated file names ----

dictionary$annot_CPRD <- gsub(".csv", ".txt", dictionary$CPRD) 
dictionary$annot_CPRD <- paste0("annot_", dictionary$annot_CPRD)
dictionary$annot_ICD <- gsub(".csv", ".txt", dictionary$ICD) 
dictionary$annot_ICD <- paste0("annot_", dictionary$annot_ICD)

write.table(dictionary, "/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
            sep = "\t", quote = F, row.names = F)

# Add read3 codes to read2 ----

merged_read_codes <- lapply(1:length(PHENOTYPES), function (i) {
  # Only do this when a file exists
  if (dictionary$CPRD[i] != "") {
    df <- read.csv(paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                          dictionary$CPRD[i]), header = T)
    # V2 code is truncated to the first 5 characters
    df$READV2_CODE <- substr(df$Readcode, 1, 5)
    # Add corresponding V3 code
    res <- merge(df, merged_v2v3[, c("READV2_CODE", "READV3_CODE")],
                 by = "READV2_CODE", all.x = T, all.y = F)
    write.table(res, 
                paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                       dictionary$annot_CPRD[i]),
                sep = "\t", quote = F, row.names = F)
  } else res <- NA
  return (res)
})

# Add ICD9 codes to ICD10 ----

merged_ICD_codes <- lapply(1:length(PHENOTYPES), function (i) {
  # Only do this when a file exists
  if (dictionary$ICD[i] != "") {
    df <- read.csv(paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                          dictionary$ICD[i]), header = T)
    # ICD code in UKB doesn't contain punctuations 
    df$orig_ICD10code <- gsub("[[:punct:]]", "", df$ICD10code)
    # Lookup file is not a perfect match - it will start with whatever ICD code
    # is in the Spiros list
    matched_res <- lapply(1:nrow(df), function (i) {
      ind <- grep(df$orig_ICD10code[i], merged_icd$ICD10)
      # When there are matches 
      if (length(ind) > 0) {
        res <- merged_icd[ind, c("ICD10", "ICD9")]
        res <- data.frame(res, df[i, ])
      } else {
        # If there are no matches just return the original Spiros code-list
        res <- df[i, ]
      }
      return (res)
    })
    res <- bind_rows(matched_res)
    res <- distinct(res)
    write.table(res,
                paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                       dictionary$annot_ICD[i]),
                sep = "\t", quote = F, row.names = F)
  } else res <- NA
  return (res)
})
