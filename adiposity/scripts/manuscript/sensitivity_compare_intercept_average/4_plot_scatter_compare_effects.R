# Author: Samvida S. Venkatesh
# Date: 23/03/23

library(tidyverse)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2211_models/GWAS/post_GWAS"
STRATA <- c("BMI_F", "BMI_M", "BMI_sex_comb",
            "Weight_F", "Weight_M", "Weight_sex_comb")

full_dat <- lapply(STRATA, function (st) {
  b0_res <- read.table(paste0(mainpath, "/gws_sig_snps/", st, "_b0_gws_sig_hits.txt"),
                       sep = "\t", header = T, stringsAsFactors = F)
  b0_res <- b0_res %>% 
    mutate(type = "b0",
           LCI = BETA - 1.96*SE,
           UCI = BETA + 1.96*SE)
  
  avg_trait_res <- read.table(paste0(mainpath, "/gws_sig_snps/", st, "_average_trait_gws_sig_hits.txt"),
                              sep = "\t", header = T, stringsAsFactors = F)
  avg_trait_res <- avg_trait_res %>% 
    mutate(type = "avg",
           LCI = BETA - 1.96*SE,
           UCI = BETA + 1.96*SE)
  
  res <- bind_rows(b0_res, avg_trait_res) %>%
    pivot_wider(names_from = type,
                values_from = c(BETA, LCI, UCI, SE, PVALUE))
  return (res)
})
names(full_dat) <- STRATA

# Plot scatter of beta_t1 vs beta_t2 ----

custom_three_diverge <- c("#D35C79","#009593", "#000000")
names(custom_three_diverge) <- c("LMM_intercept", "average_trait", "both")

plotScatter <- function (df) {
  # Decide colour
  for_plot <- df %>%
    mutate(point_type = ifelse(PVALUE_b0 <= 5E-8 & PVALUE_avg <= 5E-8, "both",
                               ifelse(PVALUE_b0 <= 5E-8, "LMM_intercept", 
                                      ifelse(PVALUE_avg <= 5E-8, "average_trait", NA)))) %>%
    filter(!is.na(point_type))
  
  # Get R2 for title
  r2_title <- cor(for_plot$BETA_b0, for_plot$BETA_avg)
  
  res_plot <- ggplot(for_plot, aes(x = BETA_b0, y = BETA_avg,
                                   fill = point_type, colour = point_type)) +
    geom_pointrange(aes(xmin = LCI_b0, xmax = UCI_b0),
                    size = 0.3, fatten = 2) +
    geom_pointrange(aes(ymin = LCI_avg, ymax = UCI_avg),
                    size = 0.3, fatten = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_abline(intercept = 0, slope = 1) +
    scale_colour_manual(values = custom_three_diverge, guide = "none") +
    scale_fill_manual(values = custom_three_diverge, guide = "none") +
    labs(title = paste0("R2 = ", round(r2_title,3))) +
    theme(axis.text = element_text(size = 6),
          axis.title = element_blank(),
          plot.title = element_text(size = 8))
  return (res_plot)
}

lapply(STRATA, function (st) {
  ggsave(paste0(mainpath, "/plots/b0_vs_avg_effect_", st, ".png"),
         plotScatter(full_dat[[st]]), 
         height = 2, width = 2, units = "in")
})

