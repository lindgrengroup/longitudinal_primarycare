# Author: Samvida S. Venkatesh
# Date: 19/02/21

library(tidyverse)

# Read relevant data ----

# UKB phenotype file 
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

# QC file from UKB
qc <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/ukb_sqc_v2.txt", header = T, 
                 na.string = c("NA", "", "."), stringsAsFactors = F)

# fam file corresponding to the QC file provided by the UKBB
fam <- read.table("/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam", 
                  header = F)

# Prepare data ----

# Add IDs to QC file
qc$eid <- fam[, 1]

# Merge QC file with phenotype file 
NPCS = 21
PCs <- paste0("PC", 1:NPCS)
pheno[, PCs] <- qc[match(pheno$eid, qc$eid), PCs]

covars <- pheno[, 
                c("eid", 
                  "f.21000.0.0", "f.21000.1.0", "f.21000.2.0",
                  "f.50.0.0", "f.50.1.0", "f.50.2.0", PCs)]
colnames(covars) <- c("eid", 
                      "ancestry.1", "ancestry.2", "ancestry.3",
                      "height.1", "height.2", "height.3", PCs)

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

covars <- cleaned 

covars$sex <- covars$Submitted.Gender

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
covars$ancestry <- apply(select(covars, starts_with("ancestry")), 1, 
                         function (x) nameAncestry(x) )

# Calculate median height
covars$height <- apply(select(covars, starts_with("height")), 1, 
                       function (x) median(as.numeric(x), na.rm = T) )

covars <- covars[, c("eid", 
                     "sex", "ancestry", "height", PCs)]

write.table(covars, "/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
            sep = "\t", quote = F, row.names = F)



