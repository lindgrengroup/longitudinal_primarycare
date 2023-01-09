# Author: Samvida S. Venkatesh
# Date: 23/10/2022

library(tidyverse)
theme_set(theme_bw())

infile_path <- "" # REDACTED
plot_dir <- "" # REDACTED

# colour palette for clusters
custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("k1", "k1_k2", "k1_k2_k3")

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
CLUSTK <- c("k1", "k1_k2", "k1_k2_k3")

# Read data ----

# Flip beta to be wrt "other_allele"

lead_snps_sumstats <- lapply(PHENOTYPES, function (pheno) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- lapply(CLUSTK, function (k) {
      
      fname <- paste0(infile_path, "/highdim_splines/standardised_outcomes/GWAS/post_GWAS/lead_snps/",
                      pheno, "_", sx, "_", k, "_all_leadsnp_assocns.txt")
      if (file.exists(fname)) {
        df <- read.table(fname,
                         sep = "\t", header = T, stringsAsFactors = F)
        
        df <- df %>% 
          mutate(phenotype = pheno,
                 sex_strata = sx,
                 cluster = k,
                 Tested_Allele = factor(Tested_Allele, levels = c("A", "C", "G", "T")),
                 Other_Allele = factor(Other_Allele, levels = c("A", "C", "G", "T")),
                 BETA_FLIP = -BETA)
      } else df <- NULL

      return (df)
    })
    res <- bind_rows(res)
    return (res)
  })
  res <- bind_rows(res_list)
  return (res)
})
lead_snps_sumstats <- bind_rows(lead_snps_sumstats)

write.table(lead_snps_sumstats,
            paste0(infile_path, "/highdim_splines/standardised_outcomes/GWAS/post_GWAS/lead_snps/all_strata_all_leadsnps_sumstats.txt"), 
            sep = "\t", quote = F, row.names = F)

# Plot results ----

lead_snps_sumstats <- read.table(paste0(infile_path, "/highdim_splines/standardised_outcomes/GWAS/post_GWAS/lead_snps/all_strata_all_leadsnps_sumstats.txt"),
                                 sep = "\t", header = T, stringsAsFactors = F)

VARIDS <- unique(lead_snps_sumstats$SNP)
VARIDS <- VARIDS[grep("^rs", VARIDS)]

plot_dat <- split(lead_snps_sumstats, f = lead_snps_sumstats$phenotype)

plot_dat <- lapply(PHENOTYPES, function (p) {
  res <- plot_dat[[p]] %>%
    mutate(BETA_FLIP = -BETA,
           sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
           cluster = factor(cluster, levels = c("k1", "k1_k2", "k1_k2_k3")),
           strata = factor(paste0(cluster, "_", sex_strata), 
                           levels = c("k1_k2_k3_sex_comb", "k1_k2_k3_M", "k1_k2_k3_F",
                                      "k1_k2_sex_comb", "k1_k2_M", "k1_k2_F",
                                      "k1_sex_comb", "k1_M", "k1_F")),
           phenotype = p,
           OR = exp(BETA_FLIP),
           uci = exp(BETA_FLIP + 1.96*SE),
           lci = exp(BETA_FLIP - 1.96*SE),
           sig_lty = factor(ifelse(PVALUE < 5E-8, "yes", "no"),
                            levels = c("yes", "no")))
  return (res)
})
names(plot_dat) <- PHENOTYPES

plotORs <- function (pheno, v) {
  sub_dat <- plot_dat[[pheno]] %>% 
    filter(SNP == v)
  
  res_plot <- ggplot(sub_dat, aes(x = OR, y = strata,
                               group = cluster)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = sex_strata, color = cluster,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.4) +
    geom_vline(xintercept = 1, linetype = 2) +
    scale_color_manual(values = custom_four_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_linetype_manual(values = c(no = 2, yes = 1)) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
 return (res_plot) 
}

lapply(PHENOTYPES, function (pheno) {
  lapply(VARIDS, function (v) {
    tiff(paste0(plot_dir, "/change_snp_effects/", pheno, "_all_clust_", v, ".tiff"),
         height = 5, width = 5, units = "cm",
         res = 300)
    print(plotORs(pheno, v))
    dev.off()
  })
})

plot_bmi_weight_k1 <- bind_rows(plot_dat) %>% 
  filter(SNP == "rs429358" & cluster == "k1") %>%
  mutate(strata = factor(paste0(phenotype, "_", sex_strata), 
                         levels = c("Weight_sex_comb", "Weight_M",
                                    "Weight_F",
                                    "BMI_sex_comb", "BMI_M", "BMI_F")))

res_plot <- ggplot(plot_bmi_weight_k1, aes(x = OR, y = strata)) +
  geom_pointrange(aes(xmin = lci, xmax = uci, 
                      shape = sex_strata, 
                      color = sex_strata,
                      linetype = sig_lty, alpha = sig_lty),
                  position = position_dodge(width = 0.7),
                  size = 0.4) +
  geom_vline(xintercept = 1, linetype = 2) +
  scale_color_manual(values = custom_three_diverge, guide = "none") +
  scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
tiff(paste0(plot_dir, "/change_snp_effects/rs429358_k1.tiff"),
     height = 5, width = 5, units = "cm",
     res = 300)
res_plot
dev.off()

