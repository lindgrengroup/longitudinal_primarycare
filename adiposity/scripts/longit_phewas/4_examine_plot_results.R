# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)
theme_set(theme_bw())

# Read data ----

res <- read.table("rs429358_all_results.txt", sep = "\t",
                  header = T, stringsAsFactors = F)

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

PHENO_LEVELS <- sort(unique(res$phenotype), decreasing = T)
STRATA_LEVELS <- paste0(rep(PHENO_LEVELS, each = 3), "_",
                        rep(c("sex_comb", "M", "F"), times = length(PHENO_LEVELS)))

PTHRESH <- 0.05/length(PHENO_LEVELS)

# Subset to full data (including discovery IDs) ----

# Full data (for supplementary plot)
plot_dat <- res %>% filter(discovery_ids_incl) %>%
  mutate(sex_strata = factor(sex_strata, levels = c("F", "M", "sex_comb")),
         phenotype = factor(phenotype, levels = PHENO_LEVELS),
         strata = paste0(phenotype, "_", sex_strata, levels = STRATA_LEVELS),
         uci = beta + se,
         lci = beta - se,
         sig_lty = factor(ifelse(pvalue < PTHRESH, "yes", "no"),
                          levels = c("yes", "no")))

MINPLOT <- min(plot_dat$lci)
MAXPLOT <- max(plot_dat$uci)

plotBetas <- function (dat) {
 
  res_plot <- ggplot(dat, aes(x = beta, y = strata,
                                  group = phenotype)) +
    geom_pointrange(aes(xmin = lci, xmax = uci,
                        shape = sex_strata, color = sex_strata,
                        linetype = sig_lty, alpha = sig_lty),
                    position = position_dodge(width = 0.7),
                    size = 0.4) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(values = custom_three_diverge, guide = "none") +
    scale_alpha_manual(values = c(no = 0.7, yes = 1)) +
    scale_x_continuous(limits = c(MINPLOT, MAXPLOT),
                       guide = guide_axis(check.overlap = TRUE)) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  return (res_plot) 
}
