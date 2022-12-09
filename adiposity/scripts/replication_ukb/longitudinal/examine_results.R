# Author: Samvida S. Venkatesh
# Date: 29/09/22

library(tidyverse)

# Read data ----

STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")
TRAITS <- c("b1", "k1", "k1_k2", "k1_k2_k3")
LEAD_SNPS <- c("rs9467663", "chr6:26076446",
               "rs11778922", "rs61955499",
               "rs12953815", "rs429358")

snp_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/"
gwas_res_dir <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/ukb_no_gp/longit_replication_GWAS/"

discovery_res <- lapply(STRATA, function (st) {
  res_list <- lapply(TRAITS, function (trt) {
    # Get file
    if (trt == "b1") {
      infile_name <- paste0(snp_dir, "2211_models/GWAS/BOLT_results/",
                            st, "_", trt, "_gws_sig_hits.txt")
    } else {
      infile_name <- paste0(snp_dir, "highdim_splines/standardised_outcomes/GWAS/BOLT_results/",
                            st, "_", trt, "_gws_sig_hits.txt")
    }
    
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
      
      res <- res %>% filter(SNP %in% LEAD_SNPS)
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
    
    res <- res %>% filter(SNP %in% LEAD_SNPS)
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
                      "SE_discovery", "SE_replication",
                      "PVALUE_discovery", "PVALUE_replication")))
    return (res)
  }
}

matched_strata_res <- lapply(STRATA, function (st) {
  res_list <- lapply(TRAITS, function (trt) {
    print(paste0(st, "_", trt))
    dres <- discovery_res[[st]][[trt]]
    rres <- replication_res[[st]][[trt]]
    if (!is.null(dres) & !is.null(rres)) {
      to_rep <- inner_join(dres, rres)
      matched_res <- flipReplicationSNPs(to_rep)
      
      # # Match direction and P < PTHRESH
      # matched_res <- matched_res %>%
      #   mutate(effect_dirn_match = 
      #            (BETA_discovery > 0 & BETA_replication > 0) | (BETA_discovery <= 0 & BETA_replication <= 0),
      #          stringent_match = effect_dirn_match & PVALUE_replication < PTHRESH,
      #          loose_match = effect_dirn_match & PVALUE_replication < 0.05,
      #          strata = paste0(st, "_", trt))
      # 
      # sink(log_file, append = T)
      # cat(paste0(st, "_", trt, ": "), "\n")
      # cat("\t", paste0("Replicated ", sum(matched_res$stringent_match), 
      #                        " of ", nrow(matched_res), " SNPs at P < ", PTHRESH), "\n")
      # cat("\t", paste0("Replicated ", sum(matched_res$loose_match), 
      #                  " of ", nrow(matched_res), " SNPs at P < 0.05"), "\n")
      # sink()
      # 
      matched_res <- matched_res %>%
        mutate(strata = paste0(st, "_", trt)) %>%
        mutate(across(all_of(c("SNP", "Tested_Allele_discovery", "Other_Allele_discovery",
                               "Tested_Allele_replication",
                               "strata")),
                      as.character)) %>%
        mutate(across(all_of(c("CHR", "POS", 
                               "BETA_discovery", "BETA_replication", 
                               "PVALUE_discovery", "PVALUE_replication")),
                      as.numeric))
      
    } else matched_res <- NULL
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

to_write <- overall_res %>% filter(PVALUE_discovery <= 5E-08)
write.table(to_write, paste0(gwas_res_dir, "lead_snps_replication_discovery_matched.txt"),
            sep = "\t", row.names = F, quote = F)

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
