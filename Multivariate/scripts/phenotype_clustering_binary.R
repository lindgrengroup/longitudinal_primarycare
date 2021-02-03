# Date: 02/02/2021
# Author: Samvida S. Venkatesh

library(tidyverse)
theme_set(theme_bw())
library(ade4)
library(pheatmap)

# Read data ----

# GP clinical annotated data
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated.rds")

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Keep demographic information
demo_info <- gp_clinical %>% distinct(eid, sex, dob, mean_UKBB_BMI)

# Read merged list of V2 and V3 codes

read_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt", 
                         sep = "\t", header = T, quote = "", fill = F,
                         comment.char = "~")

# Create binary response matrix for individual x read code ----

# Add unique code from read 2
gp_clin_new <- gp_clinical[, c("eid", "event_dt", "read_2", "read_3")]
match_ind <- match(gp_clin_new$read_2, read_codes$READV2_CODE)
gp_clin_new$unique_code_v2 <- ifelse(is.na(match_ind), NA, 
                                     read_codes$unique_code[match_ind])

# Add unique code from read 3
match_ind <- match(gp_clin_new$read_3, read_codes$READV3_CODE)
gp_clin_new$unique_code_v3 <- ifelse(is.na(match_ind), NA, 
                                     read_codes$unique_code[match_ind])

# Merge unique codes, only keeping those which are in the code-description
# file for interpretability
gp_clin_new$unique_code <- 
  ifelse(is.na(gp_clin_new$unique_code_v2) & is.na(gp_clin_new$unique_code_v3), 
         NA, ifelse(is.na(gp_clin_new$unique_code_v2), gp_clin_new$unique_code_v3,
                    gp_clin_new$unique_code_v2))

gp_clin_new <- subset(gp_clin_new, !is.na(gp_clin_new$unique_code))

# List the unique codes for each individual and pivot into a matrix

indivs <- gp_clin_new %>% group_by(eid) %>% 
  summarise(unique_codes = unique(unique_code)) %>%
  pivot_wider(id_cols = eid, names_from = unique_codes, values_from = unique_codes)

# Construct binary response matrix (code recorded or not)
response_matrix <- apply(indivs[,-1], 1, function (x) {
  ifelse(is.na(x), 0, 1)
})
colnames(response_matrix) <- as.character(indivs$eid)

# Calculate distance matrix between read codes ----

# Split by sex
SEXES <- c("F", "M")
EIDSF <- as.character(demo_info$eid[demo_info$sex == "F"])
EIDSM <- as.character(demo_info$eid[demo_info$sex == "M"])

response_matrix <- list(F = response_matrix[, EIDSF],
                        M = response_matrix[, EIDSM])

# Calculate Jaccard dissimilarity index (intersection / union metric)

jacc_dist <- lapply(response_matrix, function (mat) {
  # Remove codes for which no individual in the matrix has a measurement
  keep <- which(rowSums(mat) > 0)
  mat <- mat[keep, ]
  return (dist.binary(mat, method = 1, diag = T))
})
names(jacc_dist) <- SEXES

saveRDS(jacc_dist, "/well/lindgren/UKBIOBANK/samvida/multivariate/binary_clustering_distance_matrix.rds")



