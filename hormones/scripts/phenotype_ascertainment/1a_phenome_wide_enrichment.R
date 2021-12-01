# Author: Samvida S. Venkatesh
# Date: 29/11/2021

library(tidyverse)

# Read data ----

# Hormone "cases"
HORMONES <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/hormone_list.txt")$V1
hormone_dat <- readRDS("/well/lindgren/UKBIOBANK/samvida/full_primary_care/data/data_popn_qcd_no_longit_filter.rds")[HORMONES]
SEX_STRATA <- c("F", "M", "sex_comb")

# Annotated GP data
raw_gp_dat <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/gp_clinical_annotated.txt",
                         sep = "\t", header = T, comment.char = "$",
                         stringsAsFactors = F)

# Control ids
CONTROL_IDS <- read.table("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/data/control_ids_no_hormones_measured.txt",
                          sep = "\t", header = T, stringsAsFactors = F)

# EID x disease matrix
eid_pheno_matrix <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/eid_phenotype_matrix.txt", 
                               sep = "\t", header = T, stringsAsFactors = F)
rownames(eid_pheno_matrix) <- eid_pheno_matrix$eid
eid_pheno_matrix <- eid_pheno_matrix[, -1]

# Disease dictionary
dictionary <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F,
                         comment.char = "@", quote = "")
DISEASES <- dictionary$phenotype
colnames(eid_pheno_matrix) <- DISEASES

# Wrangle data to add sex for strata ----

hormone_dat <- lapply(hormone_dat, function (df) {
  df$sex <- raw_gp_dat$sex[match(df$eid,
                                     raw_gp_dat$eid)]
  return (df)
})

# Fisher's enrichment test ----

# Returns table of case and control frequencies, p-values for each disease

enrich_test <- function (control_ids, case_ids) {
  
  # Subset disease matrix to relevant results ----
  dis_mat <- eid_pheno_matrix[rownames(eid_pheno_matrix) %in% 
                                c(control_ids, case_ids), ]
  # Only keep diseases that have at least one patient
  # represented in the cohort
  keep_diseases <- colSums(dis_mat, na.rm = T) > 0 & 
    colSums(dis_mat, na.rm = T) < nrow(dis_mat)
  dis_mat <- dis_mat[, keep_diseases]
  DISEASES <- colnames(dis_mat)
  
  # Get disease counts and frequency in controls and cases ----
  control_counts <- colSums(dis_mat[control_ids, ], na.rm = T)
  case_counts <- colSums(dis_mat[case_ids, ], na.rm = T)
  
  counts <- data.frame(disease = DISEASES,
                       control_counts, 
                       case_counts,
                       NCONTROL = length(control_ids),
                       NCASES = length(case_ids))
  counts$control_freq <- counts$control_counts / counts$NCONTROL
  counts$case_freq <- counts$case_counts / counts$NCASES
  
  # Test for enrichment ----
  pvals <- rep(NA, length(DISEASES))
  for (d in 1:length(DISEASES)) {
    r <- counts[d, ]
    test_mat <- matrix(c(r[["case_counts"]], 
                         r[["control_counts"]],
                         r[["NCASES"]] - r[["case_counts"]],
                         r[["NCONTROL"]] - r[["control_counts"]]), 2, 2)
    pvals[d] <- fisher.test(test_mat)$p.value
  }
  
  # Report results ----
  res <- counts
  res$pval <- pvals
  res$effect_case <- ifelse(res$case_freq > res$control_freq, 
                            "over-represented", 
                            ifelse(res$case_freq < res$control_freq,
                                   "under-represented", "no change"))
  res$FDR <- p.adjust(res$pval)
  res <- res[, c("disease", "NCASES", "case_counts",
                 "pval", "FDR", "effect_case")]
  return (res)
}

# Apply enrichment test to each hormone in each strata (F, M, sex-combined) ----

disease_tables <- lapply(HORMONES, function (hr) {
  res <- lapply(SEX_STRATA, function (sx) {
    
    case_dat <- hormone_dat[[hr]]
    control_dat <- CONTROL_IDS
    
    if (sx != "sex_comb") {
      case_dat <- case_dat %>% filter(sex == sx)
      control_dat <- control_dat %>% filter(sex == sx)
    }
    
    case_ids <- as.character(unique(case_dat$eid))
    control_ids <- as.character(unique(control_dat$eid))
    
    test_result <- enrich_test(control_ids, case_ids) 
    
    # For results, report size of cohort and disease frequency in cohort
    test_result$NCOHORT <- length(case_ids) + length(control_ids)
    # Get disease counts
    dcohort <- colSums(eid_pheno_matrix[c(case_ids, control_ids), 
                                        test_result$disease],
                       na.rm = T)
    test_result$cohort_freq <- dcohort
    
    write.table(test_result, paste0("/well/lindgren/UKBIOBANK/samvida/hormone_ehr/results/",
                                hr, "/phenome_wide_enrichment_", 
                                hr, "_", sx, ".txt"),
                sep = "\t", quote = F, row.names = F)
    return ()
  })
  return ()
})
