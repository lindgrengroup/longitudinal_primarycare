# Author: Samvida S. Venkatesh
# Date: 23/10/2022

library(tidyverse)
theme_set(theme_bw())

# colour palette: 
custom_four_diverge <- c("#D35C79", "#D9AB90", "#9FBCA4", "#009593")
names(custom_four_diverge) <- c("k1", "k2", "k3", "k4")

# Read data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")
CLUSTK <- c("k1", "k2", "k3", "k4")

lead_snps_sumstats <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- lapply(CLUSTK, function (k) {
      df <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/BOLT_results/",
                              p, "_", sx, "_", k, "_all_leadsnp_assocns.txt"),
                       sep = "\t", header = T, stringsAsFactors = F)
      df$phenotype <- p
      df$sex_strata <- sx
      df$cluster <- k
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
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/highdim_splines/GWAS/BOLT_results/all_strata_all_leadsnps_sumstats.txt", 
            sep = "\t", quote = F, row.names = F)

# Plot results ----

VARIDS <- unique(lead_snps_sumstats$SNP)
VARIDS <- VARIDS[grep("^rs", VARIDS)]

plot_dat <- split(lead_snps_sumstats, f = lead_snps_sumstats$phenotype)

plot_dat <- lapply(PHENOTYPES, function (p) {
  res <- plot_dat[[p]] %>%
    mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
           cluster = factor(cluster, levels = c("k1", "k2", "k3", "k4")),
           strata = factor(paste0(cluster, "_", sex_strata), 
                           levels = c("k4_sex_comb", "k4_M", "k4_F",
                                      "k3_sex_comb", "k3_M", "k3_F",
                                      "k2_sex_comb", "k2_M", "k2_F",
                                      "k1_sex_comb", "k1_M", "k1_F")),
           OR = exp(BETA),
           uci = exp(BETA + 1.96*SE),
           lci = exp(BETA - 1.96*SE),
           sig_lty = factor(ifelse(PVALUE < 5E-8, "yes", "no"),
                            levels = c("yes", "no")))
  return (res)
})
names(plot_dat) <- PHENOTYPES

plotORs <- function (p, v) {
  sub_dat <- plot_dat[[p]] %>% 
    filter(SNP == v)
  
  MINPLOT <- min(sub_dat$lci)
  MAXPLOT <- max(sub_dat$uci)
  
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
    scale_x_continuous(limits = c(MINPLOT, MAXPLOT),
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
 return (res_plot) 
}

lapply(PHENOTYPES, function (p) {
  lapply(VARIDS, function (v) {
    tiff(paste0(p, "_all_clust_", v, ".tiff"),
         height = 4.5, width = 4.5, units = "cm",
         res = 300)
    print(plotORs(p, v))
    dev.off()
  })
})

