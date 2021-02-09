# Date: 02/02/2021
# Author: Samvida S. Venkatesh

library(tidyverse)
theme_set(theme_bw())
library(cluster)

set.seed(020221)

# Read and clean data ----

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

# Subset to keep only codes with quantitative measurement values 
gp_clinical$quant <- apply(gp_clinical[, c("value1", "value2", "value3")], 1, 
                           function (x) { any(!is.na(as.numeric(x))) } )
gp_clinical <- subset(gp_clinical, gp_clinical$quant)

# Calculate the median measurement value from all 3 value codes
gp_clinical$quant_value <- 
  apply(gp_clinical[, c("value1", "value2", "value3")], 1, 
        function (x) {
          nums <- !is.na(as.numeric(x))
          return (median(as.numeric(x[nums]))) 
          })

# Read merged list of V2 and V3 codes

read_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt", 
                         sep = "\t", header = T, quote = "", fill = F,
                         comment.char = "~")

# Add unique code from read 2
gp_clinical <- gp_clinical[, c("eid", "event_dt", "read_2", "read_3",
                               "quant_value")]
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

# IGNORING TEMPORALITY ----

## Build response matrix ----

# Calculate median measurement of each trait in each individual

indivs <- gp_clinical %>% group_by(eid, unique_code) %>% 
  summarise(median_value = median(quant_value))

response_matrix <- pivot_wider(indivs, id_cols = eid, 
                               names_from = unique_code, 
                               values_from = median_value)
EIDS <- response_matrix$eid
response_matrix <- data.matrix(response_matrix[, -1])
rownames(response_matrix) <- EIDS

## Calculate distance matrix between hormones ----

# Sex-combined

# Only keep codes that have been measured in at least 200 individuals in
# the data
keepCols <- apply(response_matrix, 2, function (x) length(which(!is.na(x))) > 200)
codenames <- ALL_CODES[keepCols]
mat <- response_matrix[, keepCols]
# Only keep individuals who have at least one of the codes measured
keepRows <- apply(mat, 1, function (x) !all(is.na(x)))
mat <- mat[keepRows, ]
# Scale each trait
mat <- t(scale(mat))
rownames(mat) <- codenames
sex_comb <- daisy(mat, metric = "gower")

# Split by sex
SEXES <- c("F", "M")
EIDSF <- demo_info$eid[demo_info$sex == "F"]

response_matrix <- list(F = response_matrix[rownames(response_matrix) %in% EIDSF, ],
                       M = response_matrix[!rownames(response_matrix) %in% EIDSF, ])

quant_dist <- lapply(SEXES, function (s) {
  
  mat <- response_matrix[[s]]
  keepCols <- apply(mat, 2, function (x) length(which(!is.na(x))) > 200)
  codenames <- ALL_CODES[keepCols]
  mat <- mat[, keepCols]
  # Only keep individuals who have at least one of the codes measured
  keepRows <- apply(mat, 1, function (x) !all(is.na(x)))
  mat <- mat[keepRows, ]
  # Scale each trait
  mat <- t(scale(mat))
  rownames(mat) <- codenames
  
  return (daisy(mat, metric = "gower"))
  
})
names(quant_dist) <- SEXES

quant_dist[["sex_comb"]] <- sex_comb

saveRDS(quant_dist, "/well/lindgren/UKBIOBANK/samvida/multivariate/quant_dist_matrix.rds")
