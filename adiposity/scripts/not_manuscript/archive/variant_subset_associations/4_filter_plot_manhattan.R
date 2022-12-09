# Author: Samvida S. Venkatesh
# Date: 22/07/22

library(tidyverse)
theme_set(theme_bw())

# Read files ----

STRATA <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/GWAS/strata_filenames.txt")$V1

ints <- lapply(STRATA, function (st) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/",
                           st, "_lmm_intercepts.txt"),
                    sep = "\t", header = T, stringsAsFactors = F,
                    comment.char = "$")
  colnames(res)[1] <- "CHROM"
  res$trait <- "intercept"
  res <- res[complete.cases(res), ]
  return (res)
})
names(ints) <- STRATA

slopes <- lapply(STRATA, function (st) {
  res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/",
                           st, "_lmm_slopes_adj_int.txt"),
                    sep = "\t", header = T, stringsAsFactors = F,
                    comment.char = "$")
  colnames(res)[1] <- "CHROM"
  res$trait <- "slope_adj_baseline"
  res <- res[complete.cases(res), ]
  return (res)
})
names(slopes) <- STRATA

# Colour palette
# light grey, light blue, teal green, amber, rose
col_palette <- c("#A4A4A4", "#A6E8F5", 
                 "#009593", "#C7B241", "#D35C79")
names(col_palette) <- c("odd_nonsig", "even_nonsig", 
                        "slope only", "intercept only", "both")

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
  int_qq <- get_qq_df(ints[[st]]$P) %>% mutate(trait = "intercept only")
  slopes_qq <- get_qq_df(slopes[[st]]$P) %>% 
    mutate(trait = "slope only")
  
  plot_qq <- bind_rows(int_qq, slopes_qq)
  
  qq_BOLT <- ggplot(plot_qq, aes(x = exp, y = obs)) +
    geom_point(aes(color = trait)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = col_palette, guide = "none") +
    labs(title = paste0("lambda intercepts = ", round(lambdaGC_ints, 3),
                        " lambda slopes = ", round(lambdaGC_slopes, 3)),
         x = "Expected -log10(P)", y = "Observed -log10(P)")
  
  ggsave(path = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/plots/",
         plot = qq_BOLT,
         filename = paste0("qq_", st, ".png"),
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
  dat <- dat %>% mutate(status = ifelse(ID %in% both_sig_ids & P < SIG_THRESH, "both",
                                       ifelse(ID %in% unique_ints & P < SIG_THRESH, "intercept only",
                                              ifelse(ID %in% unique_slopes & P < SIG_THRESH, "slope only",
                                                     ifelse(CHROM %% 2 == 0 & P >= SIG_THRESH, "even_nonsig",
                                                            ifelse(CHROM %% 2 != 0 & P >= SIG_THRESH, "odd_nonsig", NA))))))
  dat$status <- factor(dat$status)
  dat$trait <- as.factor(dat$trait)
  
  # Format data to plot
  dat <- dat %>% 
    # get chromosome length
    group_by(CHROM) %>% summarise(chr_len = max(POS)) %>%
    # get chromosome position
    mutate(tot = cumsum(as.numeric(chr_len)) - as.numeric(chr_len)) %>% select(-chr_len) %>%
    # add to original results dataset
    left_join(dat, ., by = c("CHROM" = "CHROM")) %>%
    # add cumulative position of each SNP
    arrange(CHROM, POS) %>% mutate(POS_bp = POS + tot)
  
  
  # Axis should just show chromosome number
  axisdf <- dat %>% group_by(CHROM) %>% 
    summarise(centre = (max(POS_bp) + min(POS_bp)) / 2)
  
  # Plot
  manhattan_plot <- ggplot(dat, 
                           aes(x = POS_bp, y = -log10(P)),
                           fill = status, colour = status) +
    geom_point(data = dat %>% filter(!status %in% c("both", "intercept only", "slope only")),
               aes(fill = status, colour = status), shape = 19, size = 1) +
    geom_hline(yintercept = -log10(SIG_THRESH), linetype = "dashed") +
    geom_point(data = dat %>% filter(status == "intercept only"), 
               aes(fill = status, colour = status), shape = 19, size = 1.3) +
    geom_point(data = dat %>% filter(status == "slope only"), 
               aes(fill = status, colour = status), shape = 17, size = 1.3) +
    geom_point(data = dat %>% filter(status == "both"), 
               aes(fill = status, colour = status, shape = trait), size = 1.3) +
    scale_colour_manual(values = col_palette, guide = "none") +
    scale_fill_manual(values = col_palette, guide = "none") +
    scale_shape_manual(values = c("intercept" = 19, "slope_adj_baseline" = 17),
                       guide = "none") +
    scale_x_continuous(label = axisdf$CHROM, breaks = axisdf$centre) +
    scale_y_continuous(limits = c(0, NA)) +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  
  ggsave(path = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/plots/",
         plot = manhattan_plot,
         filename = paste0("manhattan_", st, ".png"),
         device = "png", width = 10, height = 5, units = "in")
  
})

# Write results to tables ----

# Only write waist circ results for intercepts, nothing else significant anyway
to_write <- lapply(c("WC_F", "WC_M", "WC_sex_comb"), function (st) {
  SIG_THRESH <- 0.05/nrow(ints[[st]])
  to_write <- ints[[st]] %>% filter(P < SIG_THRESH) %>%
    # Align the effect sizes with the way BOLT-LMM outputs for consistency
    mutate(flip_beta = ifelse(A1 == REF, BETA, -BETA),
           strata = st, 
           chrpos = paste0(CHROM, ":", POS, ":", REF, ":", ALT),
           effect = paste0(signif(flip_beta, 3), " (", signif(SE, 3), ")")) %>%
    rename(SNP = ID) %>%
    select(all_of(c("strata", "SNP", "chrpos", "effect", "P")))
  return (to_write)
})
to_write <- bind_rows(to_write)

write.table(to_write, 
            "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/2204_models/variant_subset_associations/WC_significant_results.txt",
            sep = "\t", row.names = F, quote = F)


