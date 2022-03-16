# Author: Samvida S. Venkatesh
# Date: 19/02/21

library(tidyverse)

# Read relevant data ----

# UKB phenotype file 
pheno <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

# QC file from UKB
qc <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/QC/ukb_sqc_v2.txt", header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)

# fam file corresponding to the QC file provided by the UKBB
fam <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)

# Disease diagnoses collated from GP and UKB phenotype data
# EID x disease matrix
eid_pheno_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_phenotype_matrix.txt", 
                               sep = "\t", header = T, stringsAsFactors = F)

# Disease dictionary
dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F)
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix)[-1] <- DISEASES

# Prepare data ----

# Add IDs to QC file
qc$eid <- fam[, 1]

# Merge QC file with phenotype file 
NPCS = 21
PCs <- paste0("PC", 1:NPCS)
pheno[, PCs] <- qc[match(pheno$eid, qc$eid), PCs]

covars <- pheno[, 
                c("eid", "f.34.0.0",
                  "f.21000.0.0", "f.21000.1.0", "f.21000.2.0",
                  "f.50.0.0", "f.50.1.0", "f.50.2.0", 
                  "f.40007.0.0", "f.40007.1.0", "f.40007.2.0",
                  "f.20116.0.0", "f.20116.1.0", "f.20116.2.0",
                  PCs)]
colnames(covars) <- c("eid", "year_of_birth",
                      "ancestry.1", "ancestry.2", "ancestry.3",
                      "height.1", "height.2", "height.3", 
                      "age_at_death.1", "age_at_death.2", "age_at_death.3",
                      "smoking_status.1", "smoking_status.2", "smoking_status.3",
                      PCs)

# Merge QC file with covars
covars$Submitted.Gender <- qc$Submitted.Gender[match(covars$eid, qc$eid)]
covars$Inferred.Gender <- qc$Inferred.Gender[match(covars$eid, qc$eid)]

# Clean data ----

cleaned <- covars

# Sex mismatch
cleaned <- subset(cleaned, !is.na(cleaned$Submitted.Gender) & 
                    !is.na(cleaned$Inferred.Gender) & 
                    cleaned$Submitted.Gender == cleaned$Inferred.Gender)

# Calculate covariates ----

# Write function to name ancestry from codes
nameAncestry <- function (x) {
  ifelse(all(grepl("^1", x[!is.na(x)])), "white", 
         ifelse(all(grepl("^2", x[!is.na(x)])), "mixed",
                ifelse(all(grepl("^3", x[!is.na(x)])), "asian",
                       ifelse(all(grepl("^4", x[!is.na(x)])), "black",
                              ifelse(all(grepl("^5", x[!is.na(x)])), "chinese",
                                     ifelse(all(grepl("^6", x[!is.na(x)])), 
                                            "other",
                                            "missing or inconsistent"))))))
}

# Write function to get smoking status from codes
# 0: Never, 1: Previous, 2: Current, -3: Prefer not to answer
# Collapse to ever smoker vs never smoker (or unknown)
getSmoking <- function (x) {
  # Replace -3 with NA
  res <- x
  res[res == -3] <- NA
  # If all values are NA, assign to unknown, otherwise calculate sum of values
  # to determine: 0 - non-smoker, >0 - smoker
  smoking_status <- ifelse(all(is.na(res)), "unknown", 
                           ifelse(sum(res, na.rm = T) == 0, "non_smoker",
                                  "smoker"))
  return (smoking_status)
}

# Write function to get age at death 


# Calculate all additional covariates ----

covars <- cleaned 
covars$sex <- covars$Submitted.Gender
covars$ancestry <- apply(select(covars, starts_with("ancestry")), 1, 
                         function (x) nameAncestry(x) )
covars$smoking_status <- apply(select(covars, starts_with("smoking")), 1,
                               function (x) getSmoking(x))

# Calculate median height
covars$height <- apply(select(covars, starts_with("height")), 1, 
                       function (x) median(as.numeric(x), na.rm = T) )

# Calculate age at death (median if there are multiple values)
covars$age_at_death <- apply(select(covars, starts_with("age_at_death")),
                             1, function (x) median(as.numeric(x), na.rm = T) )

# Calculate disease status for all ICD chapters
CHAPS <- sort(unique(dictionary$ICD_chapter))

chapter_diagnosis <- lapply(CHAPS, function (chp) {
  chp_diseases <- dictionary$phenotype[dictionary$ICD_chapter == chp]
  eid_has_disease <- rowSums(eid_pheno_matrix[, chp_diseases], na.rm = T)
  return (eid_has_disease > 0)
})
eid_disease_chapters <- bind_cols(chapter_diagnosis)
colnames(eid_disease_chapters) <- paste0("ICD_chapter_", CHAPS)
eid_disease_chapters$eid <- eid_pheno_matrix$eid

covars <- merge(covars, eid_disease_chapters, by = "eid")

covars <- covars[, c("eid", "year_of_birth",
                     "sex", "ancestry", "smoking_status",
                     "height", "age_at_death",
                     paste0("ICD_chapter_", CHAPS), PCs)]

write.table(covars, "/well/lindgren/UKBIOBANK/samvida/general_resources/220131_QCd_demographic_covariates.txt",
            sep = "\t", quote = F, row.names = F)



