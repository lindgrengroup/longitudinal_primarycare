# Author: Samvida S. Venkatesh
# Date: 22/01/2021

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)
library(ggrepel)

# Read data ----

pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "f.eid"

# Remove individuals who have withdrawn consent
withdrawn <- read.table("/well/lindgren/UKBIOBANK/DATA/QC/w11867_20200820.csv", 
                        header = F)

pheno <- subset(pheno, !pheno$f.eid %in% withdrawn$V1)

case_control_ids <- readRDS("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH_FSH/case_control_ids.rds")

keep_ids <- unlist(case_control_ids)
pheno <- subset(pheno, pheno$f.eid %in% keep_ids)

# For merging ICD9 and ICD10 codes
merged_icd <- read.table("/well/lindgren/UKBIOBANK/samvida/merged_icd9_10_codes.txt",
                         sep = "\t", header = T, quote = "", fill = F)


# ICD9, ICD10, and self-reported ICD codes ----

# f.41202.x.x - main ICD10 diagnosis
# f.41204.x.x - secondary ICD10 diagnosis
# f.40001.x.x - underlying (primary) cause of death ICD10
# f.40002.x.x - contributory (secondary) cause of death ICD10
# f.41203.x.x - main ICD9 diagnosis
# f.41205.x.x - secondary ICD9 diagnosis
# f.20001.x.x - cancer code, self-reported
# f.20002.x.x - non-cancer illness code, self-reported

ICD10_cols <- c(grep("^f.41202.", colnames(pheno)),
                grep("^f.41204.", colnames(pheno)),
                grep("^f.40001.", colnames(pheno)),
                grep("^f.40002.", colnames(pheno)))
ICD10_recorded <- apply(pheno, 1, function (r) {
  codes <- r[ICD10_cols]
  codes <- merged_icd$unique_code[match(codes, merged_icd$ICD10)]
  return (codes[!is.na(codes)])
})
names(ICD10_recorded) <- pheno$f.eid

ICD9_cols <- c(grep("^f.41203.", colnames(pheno)),
               grep("^f.41205.", colnames(pheno)))
ICD9_recorded <- apply(pheno, 1, function (r) {
  codes <- r[ICD9_cols]
  codes <- merged_icd$unique_code[match(codes, merged_icd$ICD9)]
  return (codes[!is.na(codes)])
})
names(ICD9_recorded) <- pheno$f.eid

# Merge ICD9 and ICD10 codes
ICD_recorded <- mapply(c, ICD9_recorded, ICD10_recorded, SIMPLIFY = F)

cancer_cols <- grep("^f.20001.", colnames(pheno))
cancer_reported <- apply(pheno, 1, function (r) {
  codes <- r[cancer_cols]
  return (codes[!is.na(codes)])
})
names(cancer_reported) <- pheno$f.eid

non_cancer_illness_cols <- grep("^f.20002.", colnames(pheno))
non_cancer_reported <- apply(pheno, 1, function (r) {
  codes <- r[non_cancer_illness_cols]
  return (codes[!is.na(codes)])
})
names(non_cancer_reported) <- pheno$f.eid

# Calculate code counts and frequencies in cases and controls ----

# Separate by sex
sexes <- c("F", "M")
names(sexes) <- sexes
eids <- list(F = as.character(pheno$f.eid[pheno$f.31.0.0 == 0]),
             M = as.character(pheno$f.eid[pheno$f.31.0.0 == 1]))
NCASES <- lapply(sexes, function (sex) {
  length(which(case_control_ids$case %in% eids[[sex]]))
})
NCONTROLS_S1 <- lapply(sexes, function (sex) {
  length(which(case_control_ids$control_eids_s1 %in% eids[[sex]]))
})
NCONTROLS_S2 <- lapply(sexes, function (sex) {
  length(which(case_control_ids$control_eids_s2 %in% eids[[sex]]))
})
NCONTROLS_S3 <- lapply(sexes, function (sex) {
  length(which(case_control_ids$control_eids_s3 %in% eids[[sex]]))
})

counts <- lapply(list(ICD_recorded, cancer_reported,
                      non_cancer_reported), 
                 function (data) {
                   # Separate the recorded codes by sex
                   data <- list(F = data[eids[["F"]]], 
                                M = data[eids[["M"]]])
                   
                   counts <- lapply(sexes, function (sex) {
                     d <- data[[sex]]
                     cases <- unlist(d[as.character(case_control_ids$case)])
                     counts <- data.frame(table(cases))
                     colnames(counts) <- c("code", "case_count")
                     
                     control_s1 <- unlist(d[as.character(case_control_ids$control_eids_s1)])
                     counts <- merge(counts, data.frame(table(control_s1)),
                                     by.x = "code",
                                     by.y = "control_s1",
                                     all = T)
                     
                     control_s2 <- unlist(d[as.character(case_control_ids$control_eids_s2)])
                     counts <- merge(counts, data.frame(table(control_s2)),
                                     by.x = "code",
                                     by.y = "control_s2",
                                     all = T)
                     
                     control_s3 <- unlist(d[as.character(case_control_ids$control_eids_s3)])
                     counts <- merge(counts, data.frame(table(control_s3)),
                                     by.x = "code",
                                     by.y = "control_s3",
                                     all = T)
                     
                     counts[is.na(counts)] <- 0
                     colnames(counts) <- c("code", "case_count", "control_s1_count",
                                           "control_s2_count", "control_s3_count")
                     
                     # Convert counts to frequencies
                     counts$sex <- sex
                     counts$case_freq <- counts$case_count / NCASES[[sex]]
                     counts$control_s1_freq <- counts$control_s1_count / 
                       NCONTROLS_S1[[sex]]
                     counts$control_s2_freq <- counts$control_s2_count / 
                       NCONTROLS_S2[[sex]]
                     counts$control_s3_freq <- counts$control_s3_count / 
                       NCONTROLS_S3[[sex]]
                     
                     return (counts)
                   })
                   
                   return (counts)
                   
                 })
names(counts) <- c("ICD", "cancer_reported", "non_cancer_reported")

saveRDS(counts, "/well/lindgren/UKBIOBANK/samvida/hormone_ehr/LH_FSH/case_control_code_counts.rds")

# Add code descriptions ----

# ICD codes
counts$ICD <- lapply(counts$ICD, function (df) {
  # add descriptions
  df$DESC <- merged_icd$DESC[match(df$code, 
                                   merged_icd$unique_code)]
  df$ICD9 <- merged_icd$ICD9[match(df$code, 
                                   merged_icd$unique_code)]
  df$ICD10 <- merged_icd$ICD10[match(df$code, 
                                     merged_icd$unique_code)]
  return (df)
})

# Self-reported non-cancer codes (downloaded from website)
UKB_noncancer <- read.table("C:/Users/samvida/Documents/Lindgren Group/Resources/UKBIOBANK/ukb_noncancer_selfreported_codes.tsv",
                            sep = "\t", header = T, quote = "", fill = F,
                            stringsAsFactors = F)
counts$non_cancer_reported <- lapply(counts$non_cancer_reported, function (df) {
  # whitespace before and after character messes up the matching, so remove
  # whitespaces
  df$code <- trimws(df$code, which = "both")
  # add descriptions
  df$DESC <- UKB_noncancer$meaning[match(df$code, 
                                         UKB_noncancer$coding)]
  return (df)
})

# Self-reported cancer codes (downloaded from website)
UKB_cancer <- read.table("C:/Users/samvida/Documents/Lindgren Group/Resources/UKBIOBANK/ukb_cancer_selfreported_codes.tsv",
                         sep = "\t", header = T, quote = "", fill = F,
                         stringsAsFactors = F)
counts$cancer_reported <- lapply(counts$cancer_reported, function (df) {
  # whitespace before and after character messes up the matching, so remove
  # whitespaces
  df$code <- trimws(df$code, which = "both")
  # add descriptions
  df$DESC <- UKB_cancer$meaning[match(df$code, 
                                      UKB_cancer$coding)]
  return (df)
})

# Return table of top 10 diagnoses in each category of cases and controls ----

top10 <- lapply(counts, function (code_type) {
  res <- lapply(code_type, function (df) {
    keep_cases <- order(df$case_count, decreasing = T)[1:10]
    keep_control_s1 <- order(df$control_s1_count, decreasing = T)[1:10]
    keep_control_s2 <- order(df$control_s2_count, decreasing = T)[1:10]
    keep_control_s3 <- order(df$control_s3_count, decreasing = T)[1:10]
    keep <- unique(c(keep_cases, keep_control_s1, keep_control_s2,
                     keep_control_s3))
    return (df[keep, c("code", "DESC", "sex",
                       "case_freq", "control_s1_freq",
                       "control_s2_freq", "control_s3_freq")])
  })
  res <- bind_rows(res)
  return (res)
})

# Enrichment tests for codes in cases vs controls ----

enrichment_tests <- lapply(counts, function (code) {
  
  res <- lapply(code, function (df) {
    NCASES <- median(df$case_count/df$case_freq, na.rm = T)
    NCONTROLS_S1 <- median(df$control_s1_count/df$control_s1_freq, na.rm = T)
    NCONTROLS_S2 <- median(df$control_s2_count/df$control_s2_freq, na.rm = T)
    NCONTROLS_S3 <- median(df$control_s3_count/df$control_s3_freq, na.rm = T)
    
    # Fisher's tests for under- and over-representation
    
    pval_s1 <- rep(Inf, nrow(df))
    pval_s2 <- rep(Inf, nrow(df))
    pval_s3 <- rep(Inf, nrow(df))
    
    for (i in 1:nrow(df)) {
      r <- df[i, ]
      
      mat1 <- matrix(c(r[["case_count"]], 
                       r[["control_s1_count"]],
                       NCASES - r[["case_count"]],
                       NCONTROLS_S1 - r[["control_s1_count"]]), 2, 2)
      pval_s1[i] <- fisher.test(mat1)$p.value
      
      mat2 <- matrix(c(r[["case_count"]], 
                       r[["control_s2_count"]],
                       NCASES - r[["case_count"]],
                       NCONTROLS_S2 - r[["control_s2_count"]]), 2, 2)
      pval_s2[i] <- fisher.test(mat2)$p.value
      
      mat3 <- matrix(c(r[["case_count"]], 
                       r[["control_s3_count"]],
                       NCASES - r[["case_count"]],
                       NCONTROLS_S3 - r[["control_s3_count"]]), 2, 2)
      pval_s3[i] <- fisher.test(mat3)$p.value
    }
    
    df$pval_s1 <- pval_s1
    df$pval_s2 <- pval_s2
    df$pval_s3 <- pval_s3
    
    df$FDR_s1 <- p.adjust(df$pval_s1, method = "fdr")
    df$FDR_s2 <- p.adjust(df$pval_s2, method = "fdr")
    df$FDR_s3 <- p.adjust(df$pval_s3, method = "fdr")
    
    # Calculate effect direction in cases vs controls
    df$effect <- ifelse(df$case_freq > df$control_s1_freq &
                          df$case_freq > df$control_s2_freq & 
                          df$case_freq > df$control_s3_freq,
                        "over-represented", 
                        ifelse(df$case_freq < df$control_s1_freq &
                                 df$case_freq < df$control_s2_freq & 
                                 df$case_freq < df$control_s3_freq,
                               "under-represented", "inconsistent"))
    
    return (df)
  })
  
  return (res)
})

# Subset the over- and under-represented codes to print ----

enrichment_sub <- lapply(enrichment_tests, function (code_type) {
  res <- lapply(code_type, function (df) {
    # Subset codes that are nominally significantly enriched compared to all 
    # three types of controls
    # and FDR-significant in at least one
    keep <- df[which(df$pval_s1 < 0.05 & df$pval_s2 < 0.05 & df$pval_s3 < 0.05), ]
    keep <- keep[which(keep$FDR_s1 < 0.05 | keep$FDR_s2 < 0.05 | 
                         keep$FDR_s3 < 0.05), ]
    # Effect direction in cases vs controls
    #keep <- keep[keep$effect != "inconsistent", ]
    
    # Choose columns to return
    cols <- c("sex", "code", "DESC", "effect", 
              "FDR_s1", "FDR_s2", "FDR_s3",
              "case_freq", "control_s1_freq", "control_s2_freq", "control_s3_freq")
    if (any(grepl("ICD", colnames(df)))) cols <- c(cols, "ICD9", "ICD10")
    
    return (keep[, cols])
  })
  return (res)
})

# Manhattan-style plots ----

# To mirror axis, write a rough script to move the x-axis to y = 0
# Adapted from: https://stackoverflow.com/questions/39071002/moving-x-or-y-axis-together-with-tick-labels-to-the-middle-of-a-single-ggplot-n
shift_x_axis <- function(p, y = 0) {
  g <- ggplotGrob(p)
  dummy <- data.frame(y = y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, 
                                       vp = grid::viewport(y=1, 
                                                           height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
}

manhattan_plots <- lapply(enrichment_tests, function (code_type) {
  
  res <- lapply(code_type, function (df) {
    
    # Only for matched controls
    df <- df[, c("DESC", "case_freq", "control_s3_freq", "pval_s3")]
    colnames(df)[4] <- "pval"
    df$plot_pval <- -log10(df$pval)
    df$effect <- ifelse(df$case_freq > df$control_s3_freq, 
                        "over-represented in cases",
                        "under-represented in cases")
    
    # Flip "pval" for under-represented so it mirrors the y-axis
    df$plot_pval <- ifelse(df$effect == "under-represented in cases",
                            -df$plot_pval, df$plot_pval)
    
    aticks <- round(seq(min(df$plot_pval), max(df$plot_pval),
                  length.out = 5), 2)
    alabels <- abs(aticks)
    name_points_threshold <- quantile(abs(df$plot_pval), 0.999)
    
    p <- ggplot(df, aes(x = case_freq, y = plot_pval)) +
      geom_point(aes(col = effect)) +
      geom_text_repel(data = subset(df, abs(df$plot_pval) > name_points_threshold), 
                      aes(label = DESC)) +
      scale_color_manual(values = c("#A50026", "#313695"),
                         guide = F) +
      scale_y_continuous(breaks = aticks, labels = alabels) +
      labs(y = "-log10(pval)", x = "Code frequency in cases")
    
    return (shift_x_axis(p, y = 0))
  })
  return (res)
})

pdf("enriched_codes.pdf", onefile = T)
print(manhattan_plots)
dev.off()