# Date: 01/02/2021
# Author: Samvida S. Venkatesh

library(dplyr)
library(cluster)
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

# Read curated list of hormones with V2 and V3 codes

hormone_codes <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormones_v2v3_codes.txt", 
                            sep = "\t", header = T, quote = "", fill = F, 
                            comment.char = "~")

hormone_codes[hormone_codes == ""] <- NA

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

# Split by sex
SEXES <- c("F", "M")
EIDSF <- demo_info$eid[demo_info$sex == "F"]

hormone_recorded$eid <- EIDS
hormone_recorded$sex <- ifelse(hormone_recorded$eid %in% EIDSF, "F", "M")

hormone_recorded <- split(hormone_recorded, hormone_recorded$sex)
names(hormone_recorded) <- SEXES

# Calculate Jaccard dissimilarity index (intersection / union metric)

binary_dist <- lapply(hormone_recorded, function (df) {
  
  # Remove hormones (columns) for which fewer than 200 individuals in the matrix 
  # have a measurement
  eids <- df$eid
  df <- df[, HORMONES]
  keepCols <- which(colSums(df) > 200)
  
  df <- apply(df[, keepCols], 2, function (x) {as.numeric(x)} )
  
  # Remove individuals (rows) without a measurement of any of the remaining
  # hormones
  keepRows <- which(rowSums(df) > 0)
  df <- df[keepRows, ]
  
  mat <- t(df)
  
  rownames(mat) <- HORMONES[keepCols]
  colnames(mat) <- eids[keepRows]
  
  return (daisy(mat, metric = "gower", type = list(asymm = 1:ncol(mat))))
})
names(binary_dist) <- SEXES

# Cluster hormones and print ----

for (s in SEXES) {
  # Write matrix
  mat <- as.matrix(binary_dist[[s]])
  write.table(mat, paste0("binary_dist_matrix_hormones_", s, ".txt"),
              sep = "\t", quote = F, row.names = T)
  
  # Cluster (agglomerative with complete linkage)
  h <- hclust(binary_dist[[s]], method = "complete")
  
  # Print heatmaps and dendrograms
  pdf(paste0("binary_dist_matrix_hormones_", s, ".pdf"), onefile = T)
  pheatmap(mat, color = colorRampPalette(c("red", "white"))(100))
  plot(h, main = paste("Hormone clustering in", s))
  dev.off()
}