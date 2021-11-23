# Author: Samvida S. Venkatesh
# Date: 23/11/21

library(tidyverse)
theme_set(theme_bw())

set.seed(231121)

# Read files ----

STRATA <- read.table("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1
QUARTILES <- c(1:4)

assocns <- lapply(STRATA, function (st) {
  res <- lapply(QUARTILES, function (q) {
    df <- read.table(paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/PLINK_results/",
                            st, "_lmm_slopes_no_adj_b0_quartile_", q, 
                            "_assoc_results.txt"),
                     sep = "\t", header = T, stringsAsFactors = F,
                     comment.char = "$")
    colnames(df)[1] <- "CHROM"
    df$b0_quartile <- paste0("b0_q", q)
    return (df)
  })
  names(res) <- QUARTILES
  return (res)
})
names(assocns) <- STRATA

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
  
  # Within each quartile, calculate qq distribution
  all_qqs <- lapply(QUARTILES, function (q) {
    df <- get_qq_df(assocns[[st]][[q]]$P) %>% 
      mutate(b0_quartile = paste0("b0_q", q))
    return (df)
  })
  plot_qq <- bind_rows(all_qqs)
  plot_qq$b0_quartile <- factor(plot_qq$b0_quartile)
  
  qq_BOLT <- ggplot(plot_qq, aes(x = exp, y = obs)) +
    geom_point(aes(color = b0_quartile)) +
    geom_abline(intercept = 0, slope = 1) +
    labs(x = "Expected -log10(P)", y = "Observed -log10(P)")
  
  ggsave(path = paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/plots/",
                       st),
         plot = qq_BOLT,
         filename = paste0("qq_metabo_subset_b0_quartiles_", st, ".png"),
         device = "png")
})

# Manhattan plots ----

man_plots <- lapply(STRATA, function (st) {
  dat <- bind_rows(assocns[[st]])
  dat$b0_quartile <- factor(dat$b0_quartile)
  SIG_THRESH <- 0.05*4/nrow(dat) # multiply by 4 since all 4 quartiles combined

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
    geom_point(aes(color = as.factor(CHROM)), shape = 21, alpha = 0.8, size = 1) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 ), guide = F) +
    geom_point(data = subset(dat, highlight == "yes"), 
               aes(fill = b0_quartile), shape = 21, size = 2) +
    scale_fill_discrete() +
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = st) +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  ggsave(path = paste0("/well/lindgren/UKBIOBANK/samvida/adiposity/gp_only/GWAS/plots/",
                       st),
         plot = man_BOLT,
         filename = paste0("manhattan_metabo_subset_b0_quartiles_", st, ".png"),
         device = "png", width = 10, height = 5, units = "in")
  
})

