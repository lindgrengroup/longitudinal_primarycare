# Date: 02/02/2021
# Author: Samvida S. Venkatesh

library(tidyverse)
theme_set(theme_bw())
library(pheatmap)

# Read and clean data ----

# GP clinical annotated data
gp_clinical <- readRDS("/well/lindgren/UKBIOBANK/samvida/gp_clinical_annotated.rds")

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20210201.csv", 
                        header = F)
gp_clinical <- subset(gp_clinical, !gp_clinical$eid %in% withdrawn$V1)

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

# Read curated list of hormones with V2 and V3 codes

hormone_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormones_v2v3_codes.txt", 
                            sep = "\t", header = T, quote = "", fill = F,
                            comment.char = "~")

# Replace blanks with NA
hormone_codes[hormone_codes == ""] <- NA
# Assign unique code to each V2 and V3 code
hormone_codes$unique_code <- 1:nrow(hormone_codes)

# Add unique code from read 2
gp_clin_new <- gp_clinical
match_ind <- match(gp_clin_new$read_2, hormone_codes$READV2_CODE)
gp_clin_new$unique_code_v2 <- ifelse(is.na(match_ind), NA, 
                                     hormone_codes$unique_code[match_ind])

# Add unique code from read 3
match_ind <- match(gp_clin_new$read_3, hormone_codes$ctv3)
gp_clin_new$unique_code_v3 <- ifelse(is.na(match_ind), NA, 
                                     hormone_codes$unique_code[match_ind])

# Merge unique codes, only keeping those which are in the code-description
# file for interpretability
gp_clin_new$unique_code <- 
  ifelse(is.na(gp_clin_new$unique_code_v2) & is.na(gp_clin_new$unique_code_v3), 
         NA, ifelse(is.na(gp_clin_new$unique_code_v2), gp_clin_new$unique_code_v3,
                    gp_clin_new$unique_code_v2))

gp_clin_new <- subset(gp_clin_new, !is.na(gp_clin_new$unique_code))

ALL_CODES <- as.character(unique(gp_clin_new$unique_code))

# IGNORING TEMPORALITY ----

## Build response matrix ----

# Calculate median measurement of each hormone in each individual

indivs <- gp_clin_new %>% group_by(eid, unique_code) %>% 
  summarise(median_value = median(quant_value))

response_matrix <- pivot_wider(indivs, id_cols = eid, 
                               names_from = unique_code, 
                               values_from = median_value)

## Calculate distance matrix between hormones ----

# Split by sex
SEXES <- c("F", "M")
EIDSF <- demo_info$eid[demo_info$sex == "F"]

response_matrix <- list(F = response_matrix[response_matrix$eid %in% EIDSF, ],
                       M = response_matrix[!response_matrix$eid %in% EIDSF, ])

codes_sex <- lapply(response_matrix, function (df) {
  df <- df[, ALL_CODES]
  # Remove codes that have no values in any individual
  keep <- apply(df, 2, function (x) { !all(is.na(x)) })
  return (ALL_CODES[keep])
})
names(codes_sex) <- SEXES

# Calculate Euclidean distance after normalising each variable

euclidean_dist <- lapply(SEXES, function (s) {
  
  mat <- data.matrix(response_matrix[[s]][, codes_sex[[s]]])
  # Scale each hormone
  mat <- t(scale(mat))
  dist_mat <- as.matrix(dist(mat, method = "euclidean", diag = T))
  rownames(dist_mat) <- codes_sex[[s]]
  colnames(dist_mat) <- codes_sex[[s]]
  return (dist_mat)
  
})
names(euclidean_dist) <- SEXES

## Print heatmap of similarity ----

for (s in SEXES) {
  mat <- euclidean_dist[[s]]
  rownames(mat) <- hormone_codes$description[match(rownames(mat), 
                                                   hormone_codes$unique_code)]
  colnames(mat) <- hormone_codes$description[match(colnames(mat), 
                                                   hormone_codes$unique_code)]
  # Write matrix
  write.table(mat, paste0("euclidean_dist_matrix_hormones_", s, ".txt"), 
              sep = "\t", quote = F, row.names = T)
  
  # Cluster (agglomerative with complete linkage)
  h <- hclust(mat, method = "complete")
  
  # Print heatmap 
  pdf(paste0("euclidean_dist_matrix_hormones_", s, ".pdf"), onefile = T)
  pheatmap(mat, cluster_rows = F, cluster_cols = F,
           color = colorRampPalette(c("red", "white"))(100))
  plot(h, main = paste("Hormone clustering in", s))
  dev.off()
}
