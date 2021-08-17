# Author: Samvida S. Venkatesh
# Date: 31/03/21

library(tidyverse)

# Read data ----

# Read dictionary to get to phenotype lists
dictionary <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                       sep = "\t", header = T, stringsAsFactors = F,
                       quote = "")
UNIQ <- dictionary$unique_code

# Read GP clinical file
primary_care <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PRIMARY_CARE/gp_clinical.txt",
                           sep = "\t", header = T, comment.char = "$",
                           stringsAsFactors = F)
primary_care$eid <- as.character(primary_care$eid)

# Read UKB phenotype file
secondary_care <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                             header = T, sep = ",", na.string = c("NA", "", "."), 
                             stringsAsFactors = F)
colnames(secondary_care) <- gsub("X", "f.", colnames(secondary_care))
colnames(secondary_care)[1] <- "eid"

# Subset to ICD-relevant columns
# f.41202.x.x - main ICD10 diagnosis
# f.41204.x.x - secondary ICD10 diagnosis
# f.40001.x.x - underlying (primary) cause of death ICD10
# f.40002.x.x - contributory (secondary) cause of death ICD10
# f.41203.x.x - main ICD9 diagnosis
# f.41205.x.x - secondary ICD9 diagnosis

ICD9_cols <- c(grep("^f.41203.", colnames(secondary_care)),
               grep("^f.41205.", colnames(secondary_care)))
ICD10_cols <- c(grep("^f.41202.", colnames(secondary_care)),
                grep("^f.41204.", colnames(secondary_care)),
                grep("^f.40001.", colnames(secondary_care)),
                grep("^f.40002.", colnames(secondary_care)))

secondary_care <- secondary_care[, c(1, ICD9_cols, ICD10_cols)]

# Construct EID x phenotype matrix ----

ICD9_cols <- c(grep("^f.41203.", colnames(secondary_care)),
               grep("^f.41205.", colnames(secondary_care)))
ICD10_cols <- c(grep("^f.41202.", colnames(secondary_care)),
                grep("^f.41204.", colnames(secondary_care)),
                grep("^f.40001.", colnames(secondary_care)),
                grep("^f.40002.", colnames(secondary_care)))

ALL_EIDS <- unique(c(unique(primary_care$eid), secondary_care$eid))

result_matrix <- lapply(1:length(UNIQ), function (i) {
  uniq <- dictionary$unique_code[i]
  
  # Get phenotype code list for primary care (when the file exists)
  if (dictionary$CPRD[i] != "") {
    primary_care_codes <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/primary_care/",
                                            dictionary$annot_CPRD[i]), 
                                     header = T, sep = "\t", stringsAsFactors = F,
                                     quote = "", comment.char = "~")
    codes_V2 <- unique(primary_care_codes$READV2_CODE)
    codes_V2 <- codes_V2[codes_V2 != "" & !is.na(codes_V2)]
    codes_V3 <- unique(primary_care_codes$READV3_CODE)
    codes_V3 <- codes_V2[codes_V3 != "" & !is.na(codes_V3)]
    # Get IDs of individuals with these codes recorded in primary care
    match_ids_primary <- 
      unique(primary_care$eid[primary_care$read_2 %in% codes_V2 | 
                                primary_care$read_3 %in% codes_V3])
  } else match_ids_primary <- NA
  
  # Get phenotype code list for secondary care (when the file exists)
  if (dictionary$ICD[i] != "") {
    secondary_care_codes <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/secondary_care/",
                                              dictionary$annot_ICD[i]), 
                                       header = T, sep = "\t", 
                                       stringsAsFactors = F,
                                       quote = "", comment.char = "~")
    codes_ICD9 <- unique(secondary_care_codes$ICD9)
    codes_ICD9 <- codes_ICD9[codes_ICD9 != "" & !is.na(codes_ICD9)]
    codes_ICD10 <- unique(c(secondary_care_codes$ICD10, 
                            secondary_care_codes$orig_ICD10code))
    codes_ICD10 <- codes_ICD10[codes_ICD10 != "" & !is.na(codes_ICD10)]
    # Get IDs of individuals with these codes recorded in UKB
    match_ids_secondary <- apply(secondary_care, 1, function (eid_row) {
      if (any(codes_ICD9 %in% eid_row[ICD9_cols]) | 
          any(codes_ICD10 %in% eid_row[ICD10_cols])) {
        # return the ID
        res <- eid_row[1]
      } else res <- NA
      return (res)
    })
    match_ids_secondary <- unique(match_ids_secondary[!is.na(match_ids_secondary)])
  } else match_ids_secondary <- NA
  
  # Get the list of all IDs with phenotype in primary OR secondary care
  match_ids <- unique(c(match_ids_primary, match_ids_secondary))
  res <- as.numeric(ALL_EIDS %in% match_ids)
  names(res) <- ALL_EIDS
  
  return (res)
})

result_matrix <- bind_rows(result_matrix)
# transpose to get EID x disease
result_matrix <- as.data.frame(t(result_matrix))
colnames(result_matrix) <- UNIQ
result_matrix$eid <- ALL_EIDS
result_matrix <- result_matrix[, c("eid", UNIQ)]

# Save results
write.table(result_matrix, 
            "/well/lindgren/UKBIOBANK/samvida/general_resources/eid_phenotype_matrix.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
