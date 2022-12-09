# Author: Samvida S. Venkatesh
# Date: 21/05/21

library(tidyverse)
theme_set(theme_bw())

# Read files ----

SNPs_passed_QC <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/GWAS/snps_passed_QC_210520/passed_QC.txt",
                             sep = "\t", header = F, stringsAsFactors = F)

args <- commandArgs(trailingOnly = T)
STRATA <- args[1]

filename <- paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/GWAS/BOLT_results/",
                   STRATA, "_assoc.stats.gz")

if (file.exists(filename)) {
  GWAS_zip <- gzfile(filename, "rt")  
  GWAS_res <- read.table(GWAS_zip, sep = "\t", header = T, 
                         comment.char = "~", stringsAsFactors = F)
  
  # Filter GWAS results to keep only SNPs that pass QC ----
  GWAS_filtered <- subset(GWAS_res, GWAS_res$SNP %in% SNPs_passed_QC$V1)
  
  # Report metrics
  sink("/well/lindgren/UKBIOBANK/samvida/adiposity/log_files/GWAS/post_GWAS_filtering.txt", append = T)
  cat(paste0("Phenotype and Strata: ", STRATA, "\n",
             "# SNPs pre-QC: ", nrow(GWAS_res), "\n",
             "# SNPs REMAINING: ", nrow(GWAS_filtered), "\n", sep = ""))
  sink()
  
  GWAS_filtered$CHR <- as.numeric(GWAS_filtered$CHR)
  GWAS_filtered$BP <- as.numeric(GWAS_filtered$BP)
  GWAS_filtered$P_BOLT_LMM_INF <- as.numeric(GWAS_filtered$P_BOLT_LMM_INF)
  GWAS_filtered$P_LINREG <- as.numeric(GWAS_filtered$P_LINREG)
  
  saveRDS(GWAS_filtered, 
          paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/GWAS/BOLT_filtered/", 
                 STRATA, ".rds"))
  
  # Manhattan plot ----
  
  # keep only some SNPs to reduce plot size
  dat <- subset(GWAS_filtered, GWAS_filtered$P_BOLT_LMM_INF < 1E-03)
  
  dat <- dat %>% 
    # get chromosome length
    group_by(CHR) %>% summarise(chr_len = max(BP)) %>%
    # get chromosome position
    mutate(tot = cumsum(chr_len) - chr_len) %>% select(-chr_len) %>%
    # add to original results dataset
    left_join(dat, ., by = c("CHR" = "CHR")) %>%
    # add cumulative position of each SNP
    arrange(CHR, BP) %>% mutate(BP_pos = BP + tot) %>%
    # Add highlight and annotation information
    mutate(highlight = ifelse(-log10(P_BOLT_LMM_INF) > 8, "yes", "no"))
  
  # Axis should just show chromosome number
  axisdf <- dat %>% group_by(CHR) %>% 
    summarise(centre = (max(BP_pos) + min(BP_pos)) / 2)
  
  # Plot
  man_BOLT <- ggplot(dat, aes(x = BP_pos, y = -log10(P_BOLT_LMM_INF))) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
    geom_point(data = subset(dat, highlight == "yes"), 
               color = "orange", size = 2) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$centre) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = STRATA) +
    theme(legend.position = "none", 
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  pdf(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/plots/GWAS/manhattan_", 
             STRATA, ".pdf"),onefile = T)
  print(man_BOLT)
  dev.off()
}
