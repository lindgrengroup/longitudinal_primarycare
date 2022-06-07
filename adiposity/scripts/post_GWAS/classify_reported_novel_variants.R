# Author: Samvida S. Venkatesh
# Date: 17/04/2022

library(tidyverse)

# Read files ----

report_log <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/log.txt"

gp_snps <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/lmm_intercepts_sig_snps.txt",
                      sep = " ", header = F, stringsAsFactors = F)
colnames(gp_snps) <- c("SNP", "CHR", "POS")

gwas_cat_snps <- read.table("/well/lindgren/samvida/Resources/GWASCatalog/all_associations_211102.txt",
                            sep = "\t", header = T, quote = "")
gwas_cat_snps <- gwas_cat_snps %>% 
  filter(MAPPED_TRAIT == "body mass index") %>%
  select(all_of(c("CHR_ID", "CHR_POS", "SNPS"))) 

hg37_chr_pos <- read.table("/well/lindgren/samvida/Resources/grch37.all.rsid.txt",
                           sep = " ", header = T, quote = "")

# Liftover GWAS catalog chr-pos to hg19 ----

gwas_cat_snps <- gwas_cat_snps %>% 
  mutate(CHR_HG19 = hg37_chr_pos$chr[match(SNPS, hg37_chr_pos$snp)],
         POS_HG19 = hg37_chr_pos$pos[match(SNPS, hg37_chr_pos$snp)])

# Flag variants that have already been reported as BMI variants ----

sink(report_log, append = T)
cat(paste0("# Genome-wide significant SNPs: ", nrow(gp_snps), "\n"))
sink()

# Unreported variants
remove_snps <- which(gp_snps$SNP %in% gwas_cat_snps$SNPS)

gwascat_chr_pos <- paste0("chr", gwas_cat_snps$CHR_HG19, ":", gwas_cat_snps$POS)
gp_chr_pos <- paste0("chr", gp_snps$CHR, ":", gp_snps$POS)

remove_snps <- unique(c(remove_snps,
                 which(gp_chr_pos %in% gwascat_chr_pos)))

unreported_snps <- gp_snps[-remove_snps, ]

sink(report_log, append = T)
cat(paste0("\t", "Previously reported BMI SNPs: ", nrow(gp_snps) - nrow(unreported_snps), 
           "\n"))
sink()

reported_snps <- gp_snps$SNP[remove_snps]

write.table(reported_snps,
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/reported_snp_list.txt",
            sep = "\t", row.names = F, col.names = F, quote = F)

# Calculate list of conditional SNPs to assess per variant and submit GCTA-COJO ----

WINDOW_SIZE <- 500000 # 500kb window

getCondList <- function (variant_chr, variant_pos) {
  condSNPS <- gwas_cat_snps %>% filter(CHR_HG19 == variant_chr &
                                         POS_HG19 >= variant_pos - WINDOW_SIZE &
                                         POS_HG19 <= variant_pos + WINDOW_SIZE)
  return (condSNPS$SNPS)
}

for (i in 1:nrow(unreported_snps)) {
  to_write <- getCondList(unreported_snps$CHR[i],
                          unreported_snps$POS[i])
  
  # Add SNP itself to the list
  to_write <- c(unreported_snps$SNP[i], to_write)
  
  write.table(to_write,
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/reported_variants/",
                     unreported_snps$SNP[i], "_published.txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
}

submission_script <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/scripts/gcta_cojo_conditional_analysis.sh"

# Submit GCTA-COJO jobs per chromosome
for (chr in 1:22) {
  subset_chr_unreported <- unreported_snps %>% filter(CHR == chr)
  write.table(subset_chr_unreported$SNP,
              paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/GWAS/post_GWAS/BMI_sex_comb/classify_lmm_intercept_variants/tmp_variant_list_chr", 
                     chr, ".txt"),
              sep = "\t", row.names = F, col.names = F, quote = F)
  job_options <- paste(
    "-v",
    paste0(
      "CHR=\"", chr, "\""
    )
  )
  job_submission <- paste("qsub", job_options, submission_script)
  system(job_submission)
  print(job_submission)
}
