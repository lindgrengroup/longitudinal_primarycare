# Author: Samvida S. Venkatesh
# Date: 29/09/22

library(tidyverse)

# Read data ----

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")
TRAITS <- c("lmm_slopes_adj_int", "k1", "k2", "k3", "k4")
PTHRESH <- 0.05/10 ## To correct for order of magnitude of loci being tested

snp_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/"
gwas_res_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/longit_replication_GWAS/"
log_file <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/longit_replication_GWAS/snps_replicated.txt"

discovery_res <- lapply(STRATA, function (st) {
  res_list <- lapply(TRAITS, function (trt) {
    # Get file
    if (trt == "lmm_slopes_adj_int")
      infile_name <- paste0(snp_dir, "2204_models/GWAS/BOLT_results/",
                            st, "_", trt, "_gws_sig_hits.txt")
    else
      infile_name <- paste0(snp_dir, "highdim_splines/GWAS/BOLT_results/",
                            st, "_", trt, "_gws_sig_hits.txt")
    
    # Read file
    if (file.size(infile_name) == 0) 
      res <- NULL
    else {
      res <- read.table(infile_name,
                        sep = "\t", header = T, stringsAsFactors = F)
      colnames(res) <- c("SNP", "CHR", "POS", 
                         "Tested_Allele_discovery", "Other_Allele_discovery",
                         "AF_Tested_discovery", "BETA_discovery", "SE_discovery",
                         "PVALUE_discovery")
      res <- res %>%
        mutate(across(all_of(c("SNP", "Tested_Allele_discovery", "Other_Allele_discovery")),
                      as.character)) %>%
        mutate(across(all_of(c("CHR", "POS", 
                               "AF_Tested_discovery", "BETA_discovery", "SE_discovery",
                               "PVALUE_discovery")),
                      as.numeric))
    }
    return (res)
  })
  names(res_list) <- TRAITS
  return (res_list)
})
names(discovery_res) <- STRATA

replication_res <- lapply(STRATA, function (st) {
  res_list <- lapply(TRAITS, function (trt) {
    res <- read.table(paste0(gwas_res_dir, st, "/", st, "_", trt, ".txt"),
                      sep = "\t", header = T, stringsAsFactors = F, 
                      comment.char = "&")
    colnames(res) <- c("CHR", "POS", "SNP", 
                       "REF_replication", "ALT_replication", 
                       "Tested_Allele_replication",
                       "TEST_replication", "OBS_CT_replication", 
                       "BETA_replication", "SE_replication",
                       "T_STAT_replication", "PVALUE_replication")
    res <- res %>%
      mutate(across(all_of(c("SNP", 
                             "REF_replication", "ALT_replication", 
                             "Tested_Allele_replication",
                             "TEST_replication")),
                    as.character)) %>%
      mutate(across(all_of(c("CHR", "POS", 
                             "OBS_CT_replication", 
                             "BETA_replication", "SE_replication",
                             "T_STAT_replication", "PVALUE_replication")),
                    as.numeric))
    return (res)
  })
  names(res_list) <- TRAITS
  return (res_list)
}) 
names(replication_res) <- STRATA

# Match results within strata ----

flipReplicationSNPs <- function (dat) {
  res <- dat %>% 
    mutate(flag = Tested_Allele_replication != Tested_Allele_discovery,
           flip = flag & Tested_Allele_replication == Other_Allele_discovery,
           MANUAL_CHECK = flag & !flip)
  
  if (any(res$MANUAL_CHECK)) stop("Check SNPS to flip")
  else {
    res <- res %>%
      mutate(Tested_Allele_replication = ifelse(flip, Tested_Allele_discovery,
                                                Other_Allele_discovery),
             BETA_replication = ifelse(flip, -BETA_replication,
                                       BETA_replication)) %>%
      select(all_of(c("CHR", "POS", "SNP",
                      "Tested_Allele_discovery", "Other_Allele_discovery",
                      "Tested_Allele_replication",
                      "BETA_discovery", "BETA_replication",
                      "PVALUE_discovery", "PVALUE_replication")))
    return (res)
  }
}

matched_strata_res <- lapply(STRATA, function (st) {
  res_list <- lapply(TRAITS, function (trt) {
    to_rep <- inner_join(discovery_res[[st]][[trt]],
                         replication_res[[st]][[trt]])
    matched_res <- flipReplicationSNPs(to_rep)
    
    # Match direction and P < PTHRESH
    matched_res <- matched_res %>%
      mutate(effect_dirn_match = 
               (BETA_discovery > 0 & BETA_replication > 0) | (BETA_discovery <= 0 & BETA_replication <= 0),
             stringent_match = effect_dirn_match & PVALUE_replication < PTHRESH,
             loose_match = effect_dirn_match & PVALUE_replication < 0.05,
             strata = paste0(st, "_", trt))
    
    sink(log_file, append = T)
    cat(paste0(st, "_", trt, ": "), "\n")
    cat("\t", paste0("Replicated ", sum(matched_res$stringent_match), 
                           " of ", nrow(matched_res), " SNPs at P < ", PTHRESH), "\n")
    cat("\t", paste0("Replicated ", sum(matched_res$loose_match), 
                     " of ", nrow(matched_res), " SNPs at P < 0.05"), "\n")
    sink()
    
    matched_res <- matched_res %>%
      mutate(across(all_of(c("SNP", "Tested_Allele_discovery", "Other_Allele_discovery",
                             "Tested_Allele_replication",
                             "strata")),
                    as.character)) %>%
      mutate(across(all_of(c("CHR", "POS", 
                             "BETA_discovery", "BETA_replication", 
                             "PVALUE_discovery", "PVALUE_replication")),
                    as.numeric))
    
    return (matched_res)
  })
  names(res_list) <- TRAITS
  return (res_list)
})
names(matched_strata_res) <- STRATA

overall_res <- lapply(matched_strata_res, function (df_list) {
  return (bind_rows(df_list))
})
overall_res <- bind_rows(overall_res)

# Overall number of unique alleles with stringent matches
overall_stringency <- overall_res %>%
  # Get unique SNPs
  distinct(SNP, stringent_match) %>%
  # Only get best case, i.e. when there is a match
  arrange(-stringent_match) %>%
  group_by(SNP) %>% slice(1)
sum(overall_stringency$stringent_match)
mean(overall_stringency$stringent_match)

# Overall number of unique alleles with loose matches
overall_lms <- overall_res %>%
  # Get unique SNPs
  distinct(SNP, loose_match) %>%
  # Only get best case, i.e. when there is a match
  arrange(-loose_match) %>%
  group_by(SNP) %>% slice(1)
sum(overall_lms$loose_match)
mean(overall_lms$loose_match)
