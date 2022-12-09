# Author: Samvida S. Venkatesh
# Date: 23/11/21

library(tidyverse)
theme_set(theme_bw())

set.seed(231121)

# Read files ----

STRATA <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1

ints <- lapply(STRATA, function (st) {
  res <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/",
                           st, "_lmm_intercepts_assoc_results.txt"),
                    sep = "\t", header = T, stringsAsFactors = F,
                    comment.char = "$")
  colnames(res)[1] <- "CHROM"
  res$trait <- "intercept"
  return (res)
})
names(ints) <- STRATA

slopes <- lapply(STRATA, function (st) {
  res <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/",
                           st, "_lmm_slopes_adj_baseline_assoc_results.txt"),
                    sep = "\t", header = T, stringsAsFactors = F,
                    comment.char = "$")
  colnames(res)[1] <- "CHROM"
  res$trait <- "slope_adj_baseline"
  return (res)
})
names(slopes) <- STRATA

# QQ plots and lambdaGC ----

get_qq_df <- function (obs_pvals) {
  os <- obs_pvals[!is.na(obs_pvals) & obs_pvals > 0]
  os <- sort(os)
  neglog_os <- -log10(os)
  es <- (1:length(os) - 0.5)/length(os)
  neglog_es <- -log10(es)
  
  res <- data.frame(obs = neglog_os, exp = neglog_es)
  return (res)
}

qq_plots <- lapply(STRATA, function (st) {
  
  # Get genomic inflation factor 
  lambdaGC_ints <- median(qnorm(ints[[st]]$P / 2) ^ 2) / qchisq(0.5, 1)
  lambdaGC_slopes <- median(qnorm(slopes[[st]]$P / 2) ^ 2) / qchisq(0.5, 1)
  
  # Within intercepts and slopes, calculate qq distribution
  int_qq <- get_qq_df(ints[[st]]$P) %>% mutate(trait = "lmm_intercepts")
  slopes_qq <- get_qq_df(slopes[[st]]$P) %>% 
    mutate(trait = "lmm_slopes_adj_baseline")
  
  plot_qq <- bind_rows(int_qq, slopes_qq)
  
  qq_BOLT <- ggplot(plot_qq, aes(x = exp, y = obs)) +
    geom_point(aes(color = trait)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = paste0("lambda intercepts = ", round(lambdaGC_ints, 3),
                        " lambda slopes = ", round(lambdaGC_slopes, 3)),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  
  ggsave(path = paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/plots/",
                       st),
         plot = qq_BOLT,
         filename = paste0("qq_metabo_subset_", st, ".png"),
         device = "png")
})

# Manhattan plots ----

man_plots <- lapply(STRATA, function (st) {
  ints_dat <- ints[[st]] 
  slopes_dat <- slopes[[st]]
  
  SIG_THRESH <- 0.05/nrow(ints_dat)
  
  ints_sig_ids <- ints_dat$ID[ints_dat$P < SIG_THRESH]
  slopes_sig_ids <- slopes_dat$ID[slopes_dat$P < SIG_THRESH]
  both_sig_ids <- intersect(ints_sig_ids, slopes_sig_ids)
  unique_ints <- ints_sig_ids[!ints_sig_ids %in% slopes_sig_ids]
  unique_slopes <- slopes_sig_ids[!slopes_sig_ids %in% ints_sig_ids]
  
  dat <- bind_rows(ints[[st]], slopes[[st]])
  # Give points shape based on whether they cross significance
  # in both intercept-assocns and slope-assocns, or one alone
  dat <- dat %>% mutate(shape = ifelse(ID %in% both_sig_ids, "both",
                                       ifelse(ID %in% unique_ints, "intercept only",
                                              ifelse(ID %in% unique_slopes, "slope only",
                                                     "not sig"))))
  dat$shape <- factor(dat$shape)
  
  # Format data to plot
  dat <- dat %>% 
    # get chromosome length
    group_by(CHROM) %>% summarise(chr_len = max(POS)) %>%
    # get chromosome position
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
    # add to original results dataset
    left_join(dat, ., by = c("CHROM" = "CHROM")) %>%
    # add cumulative position of each SNP
    arrange(CHROM, POS) %>% mutate(POS_bp = POS + tot) %>%
    # Add highlight and annotation information
    mutate(highlight = ifelse(P < SIG_THRESH, "yes", "no"))
  
  # Axis should just show chromosome number
  axisdf <- dat %>% group_by(CHROM) %>% 
    summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)
  
  # Plot
  man_BOLT <- ggplot(dat, aes(x = POS_bp, y = -log10(P))) +
    geom_point(aes(color = as.factor(CHROM)), alpha = 0.8, size = 1.3) +
    geom_point(data = subset(dat, highlight == "yes"), 
               aes(shape = shape),
               color = "orange", size = 2) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 ), guide = F) +
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = STRATA) +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  ggsave(path = paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/plots/",
                       st),
         plot = man_BOLT,
         filename = paste0("manhattan_metabo_subset_", st, ".png"),
         device = "png", width = 10, height = 5, units = "in")
  
})

