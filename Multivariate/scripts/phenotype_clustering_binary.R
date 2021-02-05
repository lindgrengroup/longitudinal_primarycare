# Date: 02/02/2021
# Author: Samvida S. Venkatesh

library(dplyr)
library(cluster)

set.seed(020221)

# Read data ----

# GP clinical annotated data
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated.rds")

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

# Randomly subset 10,000 people for initial analyses
all_ids <- unique(gp_clinical$eid)
keep_ids <- sample(all_ids, 10000, replace = F)

gp_clinical <- subset(gp_clinical, gp_clinical$eid %in% keep_ids)

# Keep demographic information
demo_info <- gp_clinical %>% distinct(eid, sex, dob, mean_UKBB_BMI)

# Read merged list of V2 and V3 codes

read_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt", 
                         sep = "\t", header = T, quote = "", fill = F,
                         comment.char = "~")

# Create binary response matrix for individual x read code ----

# Add unique code from read 2
gp_clinical <- gp_clinical[, c("eid", "event_dt", "read_2", "read_3")]
match_ind <- match(gp_clinical$read_2, read_codes$READV2_CODE)
gp_clinical$unique_code_v2 <- ifelse(is.na(match_ind), NA, 
                                     read_codes$unique_code[match_ind])

# Add unique code from read 3
match_ind <- match(gp_clinical$read_3, read_codes$READV3_CODE)
gp_clinical$unique_code_v3 <- ifelse(is.na(match_ind), NA, 
                                     read_codes$unique_code[match_ind])

# Merge unique codes, only keeping those which are in the code-description
# file for interpretability
gp_clinical$unique_code <- 
  ifelse(is.na(gp_clinical$unique_code_v2) & is.na(gp_clinical$unique_code_v3), 
         NA, ifelse(is.na(gp_clinical$unique_code_v2), gp_clinical$unique_code_v3,
                    gp_clinical$unique_code_v2))

gp_clinical <- subset(gp_clinical, !is.na(gp_clinical$unique_code))

ALL_CODES <- unique(gp_clinical$unique_code)

EIDS <- as.character(unique(gp_clinical$eid))
eids_sex <- demo_info$sex[match(EIDS, demo_info$eid)]

EIDSF <- EIDS[eids_sex == "F"]
EIDSM <- EIDS[eids_sex == "M"]

# List the unique codes for each individual and pivot into a matrix
# of responses

indivs <- gp_clinical %>% group_by(eid) %>% 
  summarise(unique_code = as.character(unique(unique_code)))
indivs <- split(indivs[, -1], indivs$eid)

response_matrix <- lapply(indivs, function (i) {
  return (as.numeric(ALL_CODES %in% i$unique_code))
})
response_matrix <- bind_rows(response_matrix)
colnames(response_matrix) <- as.character(colnames(response_matrix))

# Calculate distance matrix between codes ----

# Sex-combined

# Only keep codes that have been measured in at least 200 individuals in
# the data
keepRows <- which(rowSums(response_matrix) > 200)
codenames <- ALL_CODES[keepRows]
mat <- response_matrix[keepRows, ]
# Only keep individuals who have at least one of the codes measured
keepCols <- which(colSums(mat) > 0)
mat <- mat[, keepCols]
rownames(mat) <- codenames
sex_comb <- daisy(mat, metric = "gower", type = list(asymm = 1:ncol(mat)))

# Split by sex
SEXES <- c("F", "M")

response_matrix <- list(F = response_matrix[, EIDSF],
                        M = response_matrix[, EIDSM])

# Calculate Gower's distance

binary_dist <- lapply(response_matrix, function (mat) {
  keepRows <- which(rowSums(mat) > 200)
  codenames <- ALL_CODES[keepRows]
  # Only keep codes that have been measured in at least 200 individuals in
  # the data
  mat <- mat[keepRows, ]
  # Only keep individuals who have at least one of the codes measured
  keepCols <- which(colSums(mat) > 0)
  mat <- mat[, keepCols]
  rownames(mat) <- codenames
  return (daisy(mat, metric = "gower", type = list(asymm = 1:ncol(mat))))
})
names(binary_dist) <- SEXES

binary_dist[["sex_comb"]] <- sex_comb

saveRDS(binary_dist, "/well/lindgren/UKBIOBANK/samvida/multivariate/binary_dist_matrix.rds")