# Author: Samvida S. Venkatesh
# Date: 22/06/2022

library(argparse)
library(tidyverse)
library(biomaRt)

parser <- ArgumentParser()
parser$add_argument("--strata", required=TRUE,
                    help = "Which strata we are assessing")
args <- parser$parse_args()

STRATA <- args$strata

filepath_main <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                        STRATA, "/classify_lmm_intercept_variants/")

# Read files ----

bm <- useMart(biomart = "ENSEMBL_MART_SNP", 
              host = "grch37.ensembl.org", path = "/biomart/martservice",
              dataset = "hsapiens_snp")

sumstats <- read.table(paste0(filepath_main, "tmp_sumstats_gcta.txt"),
                       sep = "\t", header = T, stringsAsFactors = F)
colnames(sumstats)[ncol(sumstats)] <- "SAMPLE_SIZE"
ndf_ttest <- 2*sumstats$SAMPLE_SIZE[1] - 2

potential_novel <- read.table(paste0(filepath_main, 
                                     "potential_novel_snps.txt"),
                              sep = "\t", header = F, stringsAsFactors = F)$V1

potential_refined <- read.table(paste0(filepath_main, 
                                       "potential_refined_snps.txt"),
                                sep = "\t", header = F, stringsAsFactors = F)$V1

gwas_cat_snps <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed",
                            sep = "\t", header = F, quote = "")
colnames(gwas_cat_snps) <- c("CHR", "POS0", "POS1", "SNP")

gp_snps <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/",
                             STRATA, "/finemapping/", STRATA, "_lmm_intercepts_final.lead_snps.txt"),
                      sep = " ", header = F, stringsAsFactors = F)
gp_snps <- gp_snps[, 1:3]
colnames(gp_snps) <- c("SNP", "CHR", "POS")

# Function to check if a SNP is novel ----

# Criteria for novelty: 
# 1. uncorrelated with previously reported SNPs in +/-500kb region (LD_r < 0.1)
# 2. SNP effect remains when conditioned on previously published SNPs

# RETURNS TRUE IF UNCORRELATED
isUncorr <- function (maybe_novel) {
  ld_filename <- paste0(filepath_main, "tmp_ukbb_ld/", maybe_novel,
                        "_region.ld")
  # If LD file doesn't have any data, 
  # then none of the previously reported SNPs in
  # the region are in our dataset
  if (file.size(ld_filename) != 0) {
    ld_file <- read.table(ld_filename,
                          header = F, stringsAsFactors = F)
    colnames(ld_file) <- c("CHR_A", "POS_A", "SNP_A",
                           "CHR_B", "POS_B", "SNP_B", "r2")
    return_value <- all(ld_file$r2 < 0.1)
  } else {
    return_value <- TRUE
  }
  return (return_value)
}

# RETURNS TRUE IF INDEPENDENT
ourSNPindependent <- function (our_snp) {
  file_to_check <- paste0(filepath_main, "condnl_results/", our_snp,
                          "_cond_on_pub.")
  # If the file doesn't exist, then we know SNP effect is independent 
  # of published SNPs
  # Otherwise, we have to run a check on pC (p-value after conditioning)
  if (file.exists(file_to_check)) {
    cond_res <- read.table(file_to_check, header = T, stringsAsFactors = F)
    val_return <- cond_res$pC < 0.001
  } else {
    val_return <- T
  }
  return (val_return)
}

# Functions to check if a SNP is "refined" ----

# Criteria for refinement: 
# 1. moderately correlated with previously reported SNPs in +/-500kb region (LD_r >= 0.1)
# 2. SNP effect remains when conditioned on previously published SNPs
# 3. SNP can explain all previously published correlated SNPs
# 4. SNP effect significantly stronger than the most highly correlated previously published SNP

# (1) and (2) are tested by novelty functions above

# RETURNS TRUE IF OTHER SNPs NOT INDEPENDENT
pubSNPsdependent <- function (our_snp) {
  cond_res <- read.table(paste0(filepath_main, 
                                "condnl_results/pub_cond_on_", our_snp,
                                ".cma.cojo"), header = T, stringsAsFactors = F)
  # Check if all other pCs are greater than 0.05
  return (all(cond_res$pC > 0.05))
}

# RETURNS TRUE IF OUR SNP IS STRONGER
ourSNPstronger <- function (our_snp, nearest_pub_snp) {
  # Compare effect sizes and S.E.s from GWAS results to determine
  # if our SNP has a larger effect size than the nearest published SNP
  # two-sample t-test with mean (beta) and s.e.
  our_ind <- which(sumstats$SNP == our_snp)
  pub_ind <- which(sumstats$SNP == nearest_pub_snp)
  
  mean_diff <- sumstats$BETA[our_ind] - sumstats$BETA[pub_ind]
  pooled_se <- sqrt(sumstats$SE[our_ind]^2 + sumstats$SE[pub_ind]^2)
  tstat <- mean_diff / pooled_se
  
  # Check if our SNP is significantly stronger (P < 0.05)
  return (pt(tstat, ndf_ttest) < 0.05)
}

# Functions to return nearest or most highly correlated buddy to a SNP ----

# Get buddy based on R2
getBuddyR2 <- function (lead_snp) {
  # If LD file doesn't have any entries, then no published SNP within
  # +/- 500 kb of given SNP was in our data
  ld_filename <- paste0(filepath_main, "tmp_ukbb_ld/", lead_snp,
                        "_region.ld")
  if (file.size(ld_filename) != 0) {
    ld_file <- read.table(ld_filename,
                          header = F, stringsAsFactors = F)
    colnames(ld_file) <- c("CHR_A", "POS_A", "SNP_A",
                           "CHR_B", "POS_B", "SNP_B", "r2")
    
    return_nearest <- ld_file[which.min(ld_file$r2), 
                              c("SNP_A", "CHR_A", "POS_A",
                                "SNP_B", "CHR_B", "POS_B", "r2")]
    colnames(return_nearest) <- c("SNP_og", "CHR_og", "POS_og",
                                  "SNP_buddy", "CHR_buddy", "POS_buddy", 
                                  "r2")
  } else {
    return_nearest <- data.frame(SNP_og = lead_snp,
                                 CHR_og = NA, POS_og = NA, 
                                 SNP_buddy = NA, CHR_buddy = NA, POS_buddy = NA,
                                 r2 = NA)
  }
  return (return_nearest)
}

# Get buddy based on distance
getBuddyDist <- function (lead_snp) {
  # If LD file doesn't have any entries, then no published SNP within
  # +/- 500 kb of given SNP was in our data
  ld_filename <- paste0(filepath_main, "tmp_ukbb_ld/", lead_snp,
                        "_region.ld")
  
  if (file.size(ld_filename) != 0) {
    ld_file <- read.table(ld_filename,
                          header = F, stringsAsFactors = F)
    colnames(ld_file) <- c("CHR_A", "POS_A", "SNP_A",
                           "CHR_B", "POS_B", "SNP_B", "r2")
    ld_file <- ld_file %>% 
      mutate(dist_to_novel = abs(POS_A - POS_B)) 
    
    return_nearest <- ld_file[which.min(ld_file$dist_to_novel), 
                              c("SNP_A", "CHR_A", "POS_A",
                                "SNP_B", "CHR_B", "POS_B", "dist_to_novel")]
    colnames(return_nearest) <- c("SNP_og", "CHR_og", "POS_og",
                                  "SNP_buddy", "CHR_buddy", "POS_buddy", 
                                  "dist_to_novel")
  } else {
    # Find nearest SNP from list of all obesity-associated SNPs in GWAS catalog
    chr_search <- paste0("chr", gp_snps$CHR[gp_snps$SNP == lead_snp])
    pos_test <- gp_snps$POS[gp_snps$SNP == lead_snp]
    gwascat_sub <- gwas_cat_snps %>% filter(CHR == chr_search) %>%
      mutate(dist_to_novel = abs(POS1 - pos_test))
    
    # Get nearest SNP
    ind_return <- which.min(gwascat_sub$dist_to_novel)
    return_nearest <- data.frame(SNP_og = lead_snp, 
                                 CHR_og = gp_snps$CHR[gp_snps$SNP == lead_snp],
                                 POS_og = pos_test,
                                 SNP_buddy = gwascat_sub$SNP[ind_return],
                                 CHR_buddy = gsub("chr", "", gwascat_sub$CHR[ind_return]),
                                 POS_buddy = gwascat_sub$POS1[ind_return],
                                 dist_to_novel = gwascat_sub$dist_to_novel[ind_return])
  }
  return (return_nearest)
}

# Classify potentially novel SNPs as "novel", "refined", or "reported" ----

iterate_maybe_novel <- lapply(potential_novel, function (mnsnp) {
  print(paste0("Testing SNP: ", mnsnp))
  res <- data.frame(SNP_og = mnsnp, 
                    status = "maybe novel")
  
  if (isUncorr(mnsnp) & ourSNPindependent(mnsnp)) {
    res$status <- "novel"
    res_buddy <- getBuddyDist(mnsnp)
    res <- left_join(res, res_buddy, by = "SNP_og")
  } else {
    # Is SNP stronger than the most highly correlated buddy SNP?
    res_buddy <- getBuddyR2(mnsnp)
    if (ourSNPstronger(our_snp = mnsnp, 
                       nearest_pub_snp = res_buddy$SNP_buddy)) {
      res$status <- "refined"
      res <- left_join(res, res_buddy, by = "SNP_og")
    } else {
      res$status <- "reported"
    }
  }
  return (res)
})
classified_novel <- bind_rows(iterate_maybe_novel) %>%
  filter(status != "reported")

# Classify potentially refined SNPs as "refined" or "reported" ----

iterate_maybe_refined <- lapply(potential_refined, function (mrsnp) {
  print(paste0("Testing SNP: ", mrsnp))
  res <- data.frame(SNP_og = mrsnp, 
                    status = "maybe refined")
  res_buddy <- getBuddyR2(mrsnp)
  
  if (ourSNPindependent(mrsnp) & pubSNPsdependent(mrsnp) & 
      ourSNPstronger(our_snp = mrsnp, 
                     nearest_pub_snp = res_buddy$SNP_buddy)) {
    res$status <- "refined"
    res <- left_join(res, res_buddy, by = "SNP_og")
  } else {
    res$status <- "reported"
  }
  return (res)
})
classified_refined <- bind_rows(iterate_maybe_refined) %>%
  filter(status != "reported")

# Add information from biomart to novel/refined SNPs and their buddies ----

getMAFGeneConsequenceForSNP <- function (to_annot) {
  res <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", 
                              "minor_allele_freq", 
                              "ensembl_gene_name", "consequence_type_tv"),
               filters = "snp_filter",
               values = to_annot$snpid,
               mart = bm)
  # Collapse multiple gene names and consequence types with ";"
  res <- res %>% group_by(refsnp_id) %>%
    mutate(ensembl_gene_name = paste0(unique(ensembl_gene_name), 
                                      collapse = "; "),
           consequence_type_tv = paste0(unique(consequence_type_tv), 
                                        collapse = "; ")) %>%
    distinct()
  
  return (res)
}

# Apply 

dat_to_check <- bind_rows(classified_novel,
                          classified_refined)

to_annot_og <- data.frame(snpid = dat_to_check$SNP_og,
                          chr = dat_to_check$CHR_og,
                          pos = dat_to_check$POS_og)
annot_og <- getMAFGeneConsequenceForSNP(to_annot_og)

to_annot_buddy <- data.frame(snpid = dat_to_check$SNP_buddy,
                             chr = dat_to_check$CHR_buddy,
                             pos = dat_to_check$POS_buddy)
annot_buddy <- getMAFGeneConsequenceForSNP(to_annot_buddy)

annotated_full <- dat_to_check
annotated_full <- left_join(annotated_full, annot_og,
                            by = c("SNP_og" = "refsnp_id")) %>%
  rename(MAF_og = minor_allele_freq,
         GENE_og = ensembl_gene_name,
         CONSEQUENCE_og = consequence_type_tv)
# If original data was missing chrom and pos, replace
annotated_full <- annotated_full %>%
  mutate(CHR_og = ifelse(is.na(CHR_og), chr_name, CHR_og),
         POS_og = ifelse(is.na(POS_og), chrom_start, POS_og)) %>%
  dplyr::select(-one_of(c("chr_name", "chrom_start")))

annotated_full <- left_join(annotated_full, annot_buddy,
                            by = c("SNP_buddy" = "refsnp_id")) %>%
  rename(MAF_buddy = minor_allele_freq,
         GENE_buddy = ensembl_gene_name,
         CONSEQUENCE_buddy = consequence_type_tv) 
# If original data was missing chrom and pos, replace
annotated_full <- annotated_full %>%
  mutate(CHR_buddy = ifelse(is.na(CHR_buddy), chr_name, CHR_buddy),
         POS_buddy = ifelse(is.na(POS_buddy), chrom_start, POS_buddy)) %>%
  dplyr::select(-one_of(c("chr_name", "chrom_start"))) %>%
  mutate(CHRPOS_og = paste0(CHR_og, ":", POS_og),
         CHRPOS_buddy = paste0(CHR_buddy, ":", POS_buddy))

to_write <- annotated_full %>%
  dplyr::select(all_of(c("SNP_og", "CHRPOS_og", "MAF_og", "GENE_og", "CONSEQUENCE_og",
                         "status", "dist_to_novel", "r2",
                         "SNP_buddy", "CHRPOS_buddy", "MAF_buddy", "GENE_buddy", "CONSEQUENCE_buddy")))

write.table(to_write, paste0(filepath_main,
                             "annotated_results_refined_novel_snps.txt"),
            sep = "\t", row.names = F, quote = F)

