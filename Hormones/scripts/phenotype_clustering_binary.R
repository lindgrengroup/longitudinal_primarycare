# Date: 01/02/2021
# Author: Samvida S. Venkatesh

library(tidyverse)
theme_set(theme_bw())
library(ade4)
library(pheatmap)

# Read data ----

# GP clinical annotated data
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated.rds")

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20200820.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Keep demographic information
demo_info <- gp_clinical %>% distinct(eid, sex, dob, mean_UKBB_BMI)

# Read curated list of hormones with V2 and V3 codes

hormone_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormones_v2v3_codes.txt", 
                            sep = "\t", header = T, na.strings = "")

HORMONES <- unique(hormone_codes$Group)
all_v3 <- unique(hormone_codes$ctv3)
all_v2 <- unique(hormone_codes$READV2_CODE)

# Subset gp clinical file to only those with measurement of at least
# one hormone in our list

keep <- which(gp_clinical$read_3 %in% all_v3)
keep <- c(keep, which(gp_clinical$read_2 %in% all_v2))
keep <- unique(keep)

gp_clinical <- gp_clinical[keep, ]

# Create binary response matrix for individual x hormone ----

hormone_codes <- distinct(hormone_codes, Group, ctv3, READV2_CODE)
hormone_codes <- split(hormone_codes, hormone_codes$Group)

# Create named list of V3 codes
hormones_v3 <- lapply(hormone_codes, function (df) {
  return (unique(df$ctv3))
})
# Named list of V2 codes
hormones_v2 <- lapply(hormone_codes, function (df) {
  return (unique(df$READV2_CODE))
})

EIDS <- unique(gp_clinical$eid)
hormone_recorded <- lapply(HORMONES, function (H) {
  rec_eids <- gp_clinical$eid[gp_clinical$read_3 %in% hormones_v3[[H]]]
  rec_eids <- c(rec_eids, gp_clinical$eid[gp_clinical$read_2 %in% hormones_v2[[H]]])
  rec_eids <- unique(rec_eids)
  # Return vector of 0/1 depending on whether eid is in the recorded list
  return (EIDS %in% rec_eids)
})
names(hormone_recorded) <- HORMONES

hormone_recorded <- bind_rows(hormone_recorded)

# Calculate distance matrix between hormones ----

# Jaccard index (intersection / union metric)

mat <- data.matrix(hormone_recorded)
mat <- t(mat)

jacc_dist <- dist.binary(mat, method = 1, diag = T) 

jacc_similarity <- 1 - as.matrix(jacc_dist)
write.table(jacc_similarity, "similarity_matrix_hormones.txt", sep = "\t",
            quote = F, row.names = T)

pdf("similarity_matrix_hormones.pdf")
pheatmap(jacc_similarity)
dev.off()
