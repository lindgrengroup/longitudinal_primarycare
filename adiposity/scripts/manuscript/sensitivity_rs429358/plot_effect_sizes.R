# Author: Samvida S. Venkatesh
# Date: 08/09/2022

library(tidyverse)
theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# Read data ----

plot_dir <- "" # REDACTED

PHENOS <- c("BMI", "Weight")

dat <- lapply(PHENOS, function (p) {
  res <- read.table(paste0("lme_results_", p, ".txt"),
             sep = "\t", header = T, stringsAsFactors = F)
  colnames(res) <- c("beta", "se", "tstat", "term", "adjustment",
                     "sex_strata", "phenotype")
  return (res)
})
dat <- bind_rows(dat)

# Forest plot ----

plot_dat <- dat %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         strata = factor(paste0(sex_strata, "_", adjustment),
                         levels = c("sex_comb_age_adj", "sex_comb_no_age_adj",
                                    "M_age_adj", "M_no_age_adj",
                                    "F_age_adj", "F_no_age_adj")),
         uci = beta + 1.96*se,
         lci = beta - 1.96*se)

getPlot <- function (df) {
  res_plot <- ggplot(df, aes(x = beta, y = strata)) +
    facet_wrap(~term, nrow = 1, ncol = 2, scales = "free") +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = adjustment, color = sex_strata),
                    size = 0.3,
                    position = position_dodge(width = 0.5)) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          strip.background = element_rect(size = 0.7),
          strip.text.x = element_text(size = 8))
  return (res_plot)
}

lapply(PHENOS, function (p) {
  ggsave(paste0("effect_sizes_rs429358_age_adj_", p, ".png"),
         getPlot(plot_dat %>% filter(phenotype == p)),
         units = "in", height = 3, width = 6.5)
})
