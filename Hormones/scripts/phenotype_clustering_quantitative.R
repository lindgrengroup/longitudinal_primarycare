# Date: 02/02/2021
# Author: Samvida S. Venkatesh

library(tidyverse)
library(cluster)
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

# Read V2 and V3 codes

read_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_v2v3_codes.txt",
                         sep = "\t", header = T, stringsAsFactors = F,
                         quote = "", fill = F, comment.char = "~")

# Add unique code from read 2
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

ALL_CODES <- as.character(unique(gp_clinical$unique_code))

# IGNORING TEMPORALITY ----

## Build response matrix ----

# Calculate median measurement of each hormone in each individual

indivs <- gp_clinical %>% group_by(eid, unique_code) %>% 
  summarise(median_value = median(quant_value))

response_matrix <- pivot_wider(indivs, id_cols = eid, 
                               names_from = unique_code, 
                               values_from = median_value)

## Calculate distance matrix between read codes ----

# Split by sex
SEXES <- c("F", "M")
EIDSF <- demo_info$eid[demo_info$sex == "F"]

response_matrix <- list(F = response_matrix[response_matrix$eid %in% EIDSF, ],
                       M = response_matrix[!response_matrix$eid %in% EIDSF, ])

codes_sex <- lapply(response_matrix, function (df) {
  df <- df[, ALL_CODES]
  # Remove codes that are recorded in fewer than 200 individuals
  keep <- apply(df, 2, function (x) { length(which(is.na(x))) > 200 })
  return (ALL_CODES[keep])
})
names(codes_sex) <- SEXES

# Calculate Gower distance after normalising each variable

quant_dist <- lapply(SEXES, function (s) {
  
  mat <- data.matrix(response_matrix[[s]][, codes_sex[[s]]])
  # Scale each hormone
  mat <- t(scale(mat))

  return (return(daisy(mat, metric = "gower")))
  
})
names(quant_dist) <- SEXES

## Print heatmap of distance ----

for (s in SEXES) {
  mat <- as.matrix(quant_dist[[s]])
  name_desc <- hormone_codes$description[match(codes_sex[[s]], 
                                               hormone_codes$unique_code)]
  rownames(mat) <- name_desc
  colnames(mat) <- name_desc
  # Write matrix
  write.table(mat, paste0("quant_dist_matrix_", s, ".txt"), 
              sep = "\t", quote = F, row.names = T)
  
  # Cluster (agglomerative with complete linkage)
  h <- hclust(quant_dist[[s]], method = "complete")
  
  # Print heatmap 
  pdf(paste0("quant_dist_matrix_", s, ".pdf"), onefile = T)
  pheatmap(mat, color = colorRampPalette(c("red", "white"))(100))
  plot(h, main = paste("Read code clustering in", s))
  dev.off()
}
