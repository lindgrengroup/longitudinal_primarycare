# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(tidyverse)

STRATA_ALL <- c("BMI_F", "BMI_M", "BMI_sex_comb",
                "Weight_F", "Weight_M", "Weight_sex_comb")

for (STRATA in STRATA_ALL) {
  report_log <- paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                       STRATA, "/classify_lmm_intercept_variants/log.txt")
  
  # TO RUN FOR ALL GENOME-WIDE SIGNIFICANT SNPS
  # gp_snps <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/",
  #                              STRATA, "/lmm_intercepts_sig_snps.txt"),
  #                       sep = " ", header = F, stringsAsFactors = F)
  # colnames(gp_snps) <- c("SNP", "CHR", "POS")
  
  # TO RUN FOR ONLY FINE-MAPPED SNPS
  gp_snps <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/",
                               STRATA, "/finemapping/", STRATA, "_lmm_intercepts_final.lead_snps.txt"),
                        sep = "\t", header = T, stringsAsFactors = F)
  gp_snps <- gp_snps[, 1:3]
  colnames(gp_snps) <- c("SNP", "CHR", "POS")
  
  gwas_cat_snps <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/gwascat_obesity_associations_hg19.bed",
                              sep = "\t", header = F, quote = "")
  colnames(gwas_cat_snps) <- c("CHR", "POS0", "POS1", "SNP")
  
  # Flag variants that have already been reported as obesity variants ----
  
  sink(report_log, append = T)
  cat(paste0("# Genome-wide significant SNPs: ", nrow(gp_snps), "\n"))
  sink()
  
  # Unreported variants
  remove_snps <- which(gp_snps$SNP %in% gwas_cat_snps$SNP)
  
  gwascat_chr_pos <- paste0(gwas_cat_snps$CHR, ":", gwas_cat_snps$POS1)
  gp_chr_pos <- paste0("chr", gp_snps$CHR, ":", gp_snps$POS)
  
  remove_snps <- unique(c(remove_snps,
                          which(gp_chr_pos %in% gwascat_chr_pos)))
  
  unreported_snps <- gp_snps[-remove_snps, ]
  
  sink(report_log, append = T)
  cat(paste0("\t", "Previously reported obesity SNPs: ", nrow(gp_snps) - nrow(unreported_snps), 
             "\n"))
  sink()
  
  reported_snps <- gp_snps$SNP[remove_snps]
  
  write.table(reported_snps,
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                     STRATA, "/classify_lmm_intercept_variants/reported_snp_list.txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  
  write.table(unreported_snps, 
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                     STRATA, "/classify_lmm_intercept_variants/unreported_snp_list.txt"),
              sep = "\t", row.names = F, quote = F)
  
  # Calculate list of conditional SNPs to assess per variant and submit GCTA-COJO ----
  
  WINDOW_SIZE <- 500000 # 500kb window
  
  getCondList <- function (variant_chr, variant_pos) {
    condSNPS <- gwas_cat_snps %>% filter(CHR == paste0("chr", variant_chr) &
                                           POS1 >= variant_pos - WINDOW_SIZE &
                                           POS1 <= variant_pos + WINDOW_SIZE)
    return (condSNPS$SNP)
  }
  
  for (i in 1:nrow(unreported_snps)) {
    to_write <- getCondList(unreported_snps$CHR[i],
                            unreported_snps$POS[i])
    
    # Add SNP itself to the list
    to_write <- unique(c(unreported_snps$SNP[i], to_write))
    
    write.table(to_write,
                paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                       STRATA, "/classify_lmm_intercept_variants/reported_variants/",
                       unreported_snps$SNP[i], "_published.txt"),
                sep = "\t", row.names = F, col.names = F, quote = F)
  }
  
  submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/gcta_cojo_conditional_analysis.sh"
  
  # Submit GCTA-COJO jobs per chromosome
  for (chr in 1:22) {
    subset_chr_unreported <- unreported_snps %>% filter(CHR == chr)
    
    if (nrow(subset_chr_unreported) > 0) {
      write.table(subset_chr_unreported$SNP,
                  paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/", 
                         STRATA, "/classify_lmm_intercept_variants/tmp_variant_list_chr", 
                         chr, ".txt"),
                  sep = "\t", row.names = F, col.names = F, quote = F)
      job_options <- paste(
        "-v",
        paste0(
          "STRATA=\"", STRATA, "\",",
          "CHR=\"", chr, "\""
        )
      )
      job_submission <- paste("qsub", job_options, submission_script)
      system(job_submission)
      print(job_submission)
    }
  }
}
