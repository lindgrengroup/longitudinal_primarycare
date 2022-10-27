# Author: Samvida S. Venkatesh
# Date: 21/10/22

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read data ----

PHENOTYPES <- c("BMI", "Weight")
SEX_STRATA <- c("F", "M", "sex_comb")

# WB ancestry ids
wb_ids <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eids_white_british.txt",
                     sep = "\t", header = F, stringsAsFactors = F)$V1

# Adiposity data
adipo_dat <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/indiv_qcd_data.rds")[PHENOTYPES]

# Primary care records (full)
full_gp_dat <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/gp_clinical_annotated.txt",
                          sep = "\t", header = T, comment.char = "$",
                          stringsAsFactors = F)
# Retain WB subset
full_gp_dat <- subset(full_gp_dat, full_gp_dat$eid %in% wb_ids)

# Create covariate data ----

# Summarise GP data
gp_summ <- full_gp_dat %>% 
  distinct(eid, sex, mean_UKBB_BMI) %>%
  mutate(eid = as.character(eid))

dat <- lapply(PHENOTYPES, function (p) {
  # IDs with adiposity data, subset to only primary care records
  res <- adipo_dat[[p]] %>% filter(grepl("^GP", data_provider)) 
  # Summarise to get count of number of GP visits with data and mean value across GP visits
  res <- res %>% group_by(eid) %>%
    summarise(n_visits = n(),
              mean_value = mean(value)) %>%
    mutate(eid = as.character(eid))
  # Add in individuals without GP data for this phenotype
  res <- full_join(res, gp_summ, by = "eid")
  res$n_visits[is.na(res$n_visits)] <- 0
  
  # Only keep individuals with data on sex and mean UKBB BMI
  res <- res %>% filter(!is.na(sex) & !is.na(mean_UKBB_BMI))
  return (res)
})
names(dat) <- PHENOTYPES

# Summarise for tables ----

res_tables <- lapply(PHENOTYPES, function (p) {
  full_dat <- dat[[p]] %>% 
    mutate(nvisits_group = factor(cut(n_visits, breaks = c(0, 1, 4, 10, Inf),
                                      labels = c("0", "1-3", "4-9", "10+"),
                                      right = F), 
                                  levels = c("0", "1-3", "4-9", "10+")))
  summ_res <- full_dat %>%
    group_by(nvisits_group, sex) %>%
    summarise(nindivs = n(),
              mean_mean_value = mean(mean_value),
              se_value = sd(mean_value),
              mean_mean_UKBB_BMI = mean(mean_UKBB_BMI),
              se_UKBB_BMI = sd(mean_UKBB_BMI))
  summ_res$biomarker <- p
  return (summ_res)
})
names(res_tables) <- PHENOTYPES
res_tables <- bind_rows(res_tables)

write.table(res_tables, "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/results/ascertainment_for_tables.txt",
            sep = "\t", row.names = F, quote = F)

# Re-create full adiposity (non-stratified) df for some plots
long_adiposity <- lapply(adiposity, function (x) bind_rows(x) )
names(long_adiposity) <- PHENOTYPES

# Construct demographic characteristics tables -----

demo_table <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(SS_STRATA, function (s) {
    raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "FU_n", "FUyrs", "baseline_age", "height",
          "baseline_BMI")]
  })
  df <- bind_rows(pheno)
  summary_df <- df %>% group_by(sex, ancestry) %>% 
    summarise(count = n(), 
              mean_FU_n = mean(FU_n), se_FU_n = sd(FU_n)/sqrt(count),
              mean_FUyrs = mean(FUyrs), se_FUyrs = sd(FUyrs)/sqrt(count),
              mean_bl_age = mean(baseline_age), 
              se_bl_age = sd(baseline_age)/sqrt(count),
              median_height = median(height), 
              iqr_height = paste(quantile(height, 0.25), quantile(height, 0.75),
                                 sep = ", "),
              median_bl_BMI = median(baseline_BMI), 
              iqr_bl_BMI = paste(quantile(baseline_BMI, 0.25), 
                                 quantile(baseline_BMI, 0.75),
                                 sep = ", ")) %>%
    ungroup() %>% mutate(percent = count / sum(count))
  
  adipo <- lapply(STRATA, function (s) {
    res <- adiposity[[p]][[s]][, c("eid", "age_event", "value")]
    res <- res %>% group_by(eid) %>% arrange(age_event) %>% 
      summarise(bl_value = first(value)) 
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  adipo <- bind_rows(adipo)
  adipo_summary <- adipo %>% group_by(sex, ancestry) %>% 
    summarise(median_bl_value = median(bl_value), 
              iqr_bl_value = paste(quantile(bl_value, 0.25), 
                                   quantile(bl_value, 0.75),
                                   sep = ", "))
  summary_df <- merge(summary_df, adipo_summary, by = c("sex", "ancestry"))
  
  write.table(summary_df, paste0("results/descriptive_factors/descriptive_factors_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

# Plot: ancestry and sex population counts ----

# demo_table <- lapply(PHENOTYPES, function (p) {
#   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", p, ".txt"),
#                     sep = "\t", header = T, stringsAsFactors = F)
# })
# names(demo_table) <- PHENOTYPES

popn_ancestry_treemap <- lapply(PHENOTYPES, function (p) {
  df <- demo_table[[p]] %>% group_by(ancestry) %>% summarise(size = sum(count))
  tiff(filename = paste0("plots/descriptive_factors/popn_ancestry_", p, ".tiff"),
       height = 10, width = 10, units = "cm",
       res = 600)
  treemap(df, 
          index = "ancestry", vSize = "size", type = "index",
          palette = "Dark2", border.col = c("black"), border.lwds = 1,
          title = "", fontsize.labels = 0)
  dev.off()
  return ()
})

popn_sex_pies <- lapply(PHENOTYPES, function (p) {
  df <- demo_table[[p]]
  for (a in unique(df$ancestry)) {
    res <- ggplot(subset(df, df$ancestry == a), 
                  aes(x = "", y = count, fill = sex)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar("y", start = 0) +
      scale_color_brewer(palette = "Dark2") +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.title.y = element_blank(), panel.border = element_blank(), 
            panel.grid = element_blank(), axis.ticks = element_blank(),
            legend.position = "none")
    tiff(filename = paste0("plots/descriptive_factors/popn_sex_", p, "_", 
                           a, ".tiff"),
         height = 10, width = 10, units = "cm",
         res = 600)
    print(res)
    dev.off()
  }
  return ()
})

# Box and violin plots of distributions ----

plot_distributions <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(SS_STRATA, function (s) {
    raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "FU_n", "FUyrs", 
                             "baseline_age", "height", 
                             "baseline_BMI", "baseline_trait")]
  })
  df <- bind_rows(pheno)
  
  # Number of follow-up measures
  FU_n <- ggplot(df, aes(x = ancestry, y = FU_n, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "# repeat measures") +
    theme(legend.position = "none")
  
  # Number of years of follow up
  FUyrs <- ggplot(df, aes(x = ancestry, y = FUyrs, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "follow-up length (years)") +
    theme(legend.position = "none")
  
  # Baseline age
  bl_age <- ggplot(df, aes(x = ancestry, y = baseline_age)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(20, 80), breaks = seq(20, 80, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline age (years)") +
    theme(legend.position = "none")
  
  # Adult height
  med_height <- ggplot(df, aes(x = ancestry, y = height)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(110, 210), breaks = seq(110, 210, by = 10)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "median adult height (cm)") +
    theme(legend.position = "none")
  
  # Baseline BMI
  bl_BMI <- ggplot(df, aes(x = ancestry, y = baseline_BMI)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_y_continuous(limits = c(10, 65), breaks = seq(10, 65, by = 5)) +
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline BMI (kg/m2)") +
    theme(legend.position = "none")
  
  # Baseline adiposity trait
  bl_trait <- ggplot(df, aes(x = ancestry, y = baseline_trait)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "baseline adiposity trait") +
    theme(legend.position = "none")
  
  pdf(paste0("plots/descriptive_factors/covariate_distributions_", p, ".pdf"),
      onefile = T)
  print(FU_n)
  print(FUyrs)
  print(bl_age)
  print(med_height)
  print(bl_BMI)
  print(bl_trait)
  dev.off()
  return ()
})

# Summary plot of baseline age, FU length, FU measures -----

# demo_table <- lapply(PHENOTYPES, function (p) {
#   res <- read.table(paste0("results/descriptive_factors/descriptive_factors_", 
#                            p, ".txt"),
#                     sep = "\t", header = T, stringsAsFactors = F)
# })
# names(demo_table) <- PHENOTYPES

demo_table <- lapply(PHENOTYPES, function (p) {
  res <- demo_table[[p]]
  res$adipo_trait <- p
  return (res)
})
demo_table <- bind_rows(demo_table)
# min-max normalise line thickness from 1:5
demo_table$line_thick <- 1 + 
  4*(demo_table$mean_FU_n - min(demo_table$mean_FU_n)) / 
  (max(demo_table$mean_FU_n) - min(demo_table$mean_FU_n))

pdf("plots/descriptive_factors/summary_FU.pdf")
ggplot(demo_table) +
  facet_wrap(~ancestry+sex, ncol = 2) +
  geom_segment(aes(x = mean_bl_age, xend = mean_bl_age + mean_FUyrs,
                   y = adipo_trait, yend = adipo_trait, 
                   size = mean_FU_n,
                   color = ancestry)) +
  scale_y_discrete(limits = c("WHR", "waist", "weight", "BMI"),
                   labels = c("WHR", "WC", "weight", "BMI")) +
  scale_x_continuous(limits = c(40, 65), breaks = seq(40, 65, by = 5)) +
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "Age (years)", y = "Adiposity trait") +
  theme(legend.position = "none", 
        strip.background = element_blank(),
        strip.text.x = element_blank())
dev.off()

# Plot binned trajectories ----

# Calculate mean and SE in each 5-year interval in each strata

age_bin_cuts <- seq(20, 80, by = 5)

binned_adipo <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    df <- adiposity[[p]][[s]]
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    res <- df %>% group_by(age_bin) %>% summarise(count = n(),
                                                  mean_value = mean(value),
                                                  se_value = sd(value)/sqrt(count))
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  res <- bind_rows(res)
  return (res)
})
names(binned_adipo) <- PHENOTYPES

traj_plots <- lapply(PHENOTYPES, function (p) {
  res <- ggplot(binned_adipo[[p]], aes(x = age_bin, y = mean_value,
                                       group = ancestry,
                                       color = ancestry, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_point() +
    geom_path() +
    geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                alpha = 0.2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
         title = p)
  return (res)
})
names(traj_plots) <- PHENOTYPES

pdf("plots/descriptive_factors/binned_mean_trajectories.pdf", onefile = T)
print(traj_plots)
dev.off()

# Plot sample trajectories ----

sample_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    df <- adiposity[[p]][[s]]
    ids <- sample(df$eid, 10, replace = F)
    
    plot_df <- long_adiposity[[p]]
    plot_df <- subset(plot_df, plot_df$eid %in% ids)
    
    p <- ggplot(plot_df, aes(x = age_event, y = value, group = eid)) +
      geom_line() +
      geom_point() +
      labs(x = "age (years)", y = "Adiposity trait value", 
           title = s)
    return (p)
  })
  names(res) <- SS_STRATA
  pdf(paste0("plots/descriptive_factors/sample_trajectories_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(sample_trajectories) <- PHENOTYPES

# Plot change in trait vs mean trait (sample) ----

change_mean_dfs <- lapply(PHENOTYPES, function (pheno) {
  res <- long_adiposity[[pheno]] %>% group_by(eid) %>%
    arrange(age_event, .by_group = T) %>% 
    mutate(mean_value_t1t2 = (value + lag(value))/2,
           value_change = value - lag(value),
           age_change = age_event - lag(age_event),
           value_change_t1t2 = value_change / age_change)
})

sample_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    df <- adiposity[[p]][[s]]
    ids <- sample(df$eid, 10, replace = F)
    
    plot_df <- change_mean_dfs[[p]]
    plot_df <- subset(plot_df, plot_df$eid %in% ids)
    
    p <- ggplot(plot_df, aes(x = mean_value_t1t2, 
                             y = value_change_t1t2, group = eid)) +
      geom_line() +
      geom_point() +
      labs(x = "Mean adiposity trait value", y = "Change in adiposity / time", 
           title = s)
    return (p)
  })
  names(res) <- SS_STRATA
  pdf(paste0("plots/descriptive_factors/sample_change_v_mean_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(sample_trajectories) <- PHENOTYPES

# Plot change in trait vs mean trait (binned) ----

binned_adipo <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    # Calculate mean and SE in each trait bin in each strata
    df <- change_mean_dfs[[p]]
    trait_bin_cuts <- seq(min(df$mean_value_t1t2, na.rm = T), 
                          max(df$mean_value_t1t2, na.rm = T), 
                          length.out = 11)
    
    df$trait_bin <- cut(df$mean_value_t1t2, trait_bin_cuts, include.lowest = T)
    res <- df %>% group_by(trait_bin) %>% 
      summarise(count = n(),
                mean_change_value = mean(value_change_t1t2, na.rm = T),
                se_value = sd(value_change_t1t2)/sqrt(count))
    res$ancestry <- strsplit(s, "_")[[1]][1]
    res$sex <- strsplit(s, "_")[[1]][2]
    return (res)
  })
  res <- bind_rows(res)
  return (res)
})
names(binned_adipo) <- PHENOTYPES

traj_plots <- lapply(PHENOTYPES, function (p) {
  
  res <- ggplot(binned_adipo[[p]], aes(x = age_bin, y = mean_value,
                                       group = ancestry,
                                       color = ancestry, fill = ancestry)) +
    facet_wrap(~sex, nrow = 2) +
    geom_point() +
    geom_path() +
    geom_ribbon(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                alpha = 0.2) +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
         title = p)
  return (res)
})
names(traj_plots) <- PHENOTYPES

pdf("plots/descriptive_factors/binned_mean_trajectories.pdf", onefile = T)
print(traj_plots)
dev.off()

# Construct raw slope summary table ----

rs_table <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(SS_STRATA, function (s) {
    raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "raw_slope")]
  })
  df <- bind_rows(pheno)
  summary_df <- df %>% group_by(sex, ancestry) %>% 
    summarise(median_slope = median(raw_slope), 
              iqr_slope = paste(quantile(raw_slope, 0.25), 
                                quantile(raw_slope, 0.75),
                                sep = ", "))
  write.table(summary_df, paste0("results/raw_slopes/raw_slope_summaries_", p, ".txt"),
              sep = "\t", quote = F, row.names = F)
  return (summary_df)
})

# Plot violin distribution of raw slopes ----

rs_distributions <- lapply(PHENOTYPES, function (p) {
  pheno <- lapply(SS_STRATA, function (s) {
    raw_slopes[[p]][[s]][, c("eid", "sex", "ancestry", "raw_slope")]
  })
  df <- bind_rows(pheno)
  res <- ggplot(df, aes(x = ancestry, y = raw_slope)) +
    facet_wrap(~sex, nrow = 2) +
    geom_violin(aes(fill = ancestry), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_x_discrete(limits = 
                       c("white", "asian", "other", "black", "mixed", "chinese")) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "ancestry", y = "raw slope", title = p) +
    theme(legend.position = "none")
  return (res)
})
pdf(paste0("plots/raw_slopes/raw_slopes_distributions.pdf"),
    onefile = T)
print(rs_distributions)
dev.off()

# Plot trajectories of sample individuals at the tails of raw slope distribution ----

rs_tail_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    rs_df <- raw_slopes[[p]][[s]]
    top_ids <- rs_df$eid[rs_df$raw_slope > quantile(rs_df$raw_slope, 0.95)]
    top_ids <- sample(top_ids, min(length(top_ids), 5), replace = F)
    bottom_ids <- rs_df$eid[rs_df$raw_slope < quantile(rs_df$raw_slope, 0.05)]
    bottom_ids <- sample(bottom_ids, min(length(bottom_ids), 5), replace = F)
    
    plot_df <- long_adiposity[[p]]
    plot_df <- subset(plot_df, 
                      plot_df$eid %in% top_ids | plot_df$eid %in% bottom_ids)
    plot_df$tail <- ifelse(plot_df$eid %in% top_ids, "top", "bottom")
    
    p <- ggplot(plot_df, aes(x = age_event, y = value, 
                             group = eid, color = tail)) +
      geom_line() +
      geom_point() +
      scale_color_manual(values = c("top" = "#984EA3", "bottom" = "#E41A1C"), 
                         guide = F) +
      labs(x = "age (years)", y = "Adiposity trait value", 
           title = s)
    return (p)
  })
  names(res) <- SS_STRATA
  pdf(paste0("plots/raw_slopes/tail_trajectories_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(rs_tail_trajectories) <- PHENOTYPES

# Plot mean (binned) trajectories by raw slope quartile ----

rs_quartile_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(SS_STRATA, function (s) {
    rs <- raw_slopes[[p]][[s]]
    rs$q <- cut(rs$raw_slope, quantile(rs$raw_slope), include.lowest = T,
                labels = paste0("q", 1:4))
    df <- long_adiposity[[p]]
    df <- subset(df, df$eid %in% rs$eid)
    df$q <- rs$q[match(df$eid, rs$eid)]
    
    # Calculate mean and SE in each 5-year interval within each quartile
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    plot_df <- df %>% group_by(q, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    # Plot 
    p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                             group = q, color = q, fill = q)) +
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - se_value, 
                      ymax = mean_value + se_value),
                  alpha = 0.2) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
           title = s)
    return (p)
  })
  names(res) <- SS_STRATA
  pdf(paste0("plots/raw_slopes/mean_trajectories_quartile_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(rs_quartile_trajectories) <- PHENOTYPES

# Plot correlations between raw slope and model covariates ----

rs_correlations <- lapply(PHENOTYPES, function (p) {
  
  res <- lapply(SC_STRATA, function (s) {
    df <- raw_slopes[[p]][[s]][, c("eid", "raw_slope", 
                             "sex", "ancestry", 
                             "FU_n", "FUyrs", 
                             "baseline_age", "height", 
                             "baseline_BMI", "baseline_trait")]
    # Sample 1000 IDs to plot because otherwise too many to interpret
    SAMPLE_IDS <- sample(df$eid, min(1000, nrow(df)), replace = F)
    # But calculate correlations from all
    plot_df <- pivot_longer(df, cols = c(FU_n, FUyrs, baseline_age, height,
                                         baseline_BMI, baseline_trait),
                            names_to = "covariate", values_to = "covariate_value")
    point_plots <- subset(plot_df, plot_df$eid %in% SAMPLE_IDS)
    
    # Plot correlations 
    res <- ggplot(plot_df, aes(x = covariate_value, y = raw_slope,
                               color = sex, fill = sex)) +
      facet_wrap(~covariate, ncol = 2, scales = "free_x") + 
      geom_point(data = point_plots, 
                 aes(x = covariate_value, y = raw_slope, color = sex),
                 alpha = 0.5) +
      geom_smooth(method = "lm") +
      labs(x = "Covariate value", y = "Raw slope", 
           title = paste(unique(df$ancestry), "ancestry"))
    
    return (res)
  })
  names(res) <- SC_STRATA
  pdf(paste0("plots/raw_slopes/covariate_correlations_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(rs_correlations) <- PHENOTYPES


# DEBUGGING: # Plot characteristics of q1 vs q4 by raw slope, white ----

rs_white_descriptive <- lapply(PHENOTYPES, function (p) {
  df <- lapply(c("white_F", "white_M"), function (s) {
    rs <- raw_slopes[[p]][[s]]
    rs$q <- cut(rs$raw_slope, quantile(rs$raw_slope), include.lowest = T,
                labels = paste0("q", 1:4))
    rs <- subset(rs, rs$q %in% c("q1", "q4"))
    return (rs)
  })
  df <- bind_rows(df)
  
  # Number of follow-up measures
  FU_n <- ggplot(df, aes(x = q, y = FU_n, fill = q)) +
    facet_wrap(~sex, ncol = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
    labs(x = "raw slope quartile", y = "# repeat measures") +
    theme(legend.position = "none")
  
  # Number of years of follow up
  FUyrs <- ggplot(df, aes(x = q, y = FUyrs, fill = q)) +
    facet_wrap(~sex, ncol = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
    labs(x = "raw slope quartile", y = "# follow-up years") +
    theme(legend.position = "none")
  
  # Baseline age
  bl_age <- ggplot(df, aes(x = q, y = baseline_age, fill = q)) +
    facet_wrap(~sex, ncol = 2) +
    geom_boxplot(position = position_dodge(1)) +
    scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
    labs(x = "raw slope quartile", y = "baseline age") +
    theme(legend.position = "none")
  
  # Baseline BMI
  bl_BMI <- ggplot(df, aes(x = q, y = baseline_BMI)) +
    facet_wrap(~sex, ncol = 2) +
    geom_violin(aes(fill = q), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
    labs(x = "raw slope quartile", y = "baseline BMI (kg/m2)") +
    theme(legend.position = "none")
  
  # Baseline adiposity trait
  bl_trait <- ggplot(df, aes(x = q, y = baseline_trait)) +
    facet_wrap(~sex, ncol = 2) +
    geom_violin(aes(fill = q), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_manual(values = c("q4" = "#984EA3", "q1" = "#E41A1C")) + 
    labs(x = "raw slope quartile", y = "baseline trait") +
    theme(legend.position = "none")
  
  pdf(paste0("plots/raw_slopes/white_q1_v_q4_covariate_distributions_", p, ".pdf"),
      onefile = T)
  print(FU_n)
  print(FUyrs)
  print(bl_age)
  print(bl_BMI)
  print(bl_trait)
  dev.off()
  return ()
})

# Plot trajectories of sample individuals at the tails of adj slope distribution ----

as_tail_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    rs_df <- adj_slopes[[p]][[s]]
    top_ids <- rs_df$eid[rs_df$residual > quantile(rs_df$residual, 0.95)]
    top_ids <- sample(top_ids, min(length(top_ids), 5), replace = F)
    bottom_ids <- rs_df$eid[rs_df$residual < quantile(rs_df$residual, 0.05)]
    bottom_ids <- sample(bottom_ids, min(length(bottom_ids), 5), replace = F)
    
    plot_df <- long_adiposity[[p]]
    plot_df <- subset(plot_df, 
                      plot_df$eid %in% top_ids | plot_df$eid %in% bottom_ids)
    plot_df$tail <- ifelse(plot_df$eid %in% top_ids, "top", "bottom")
    
    p <- ggplot(plot_df, aes(x = age_event, y = value, 
                             group = eid, color = tail)) +
      geom_line() +
      geom_point() +
      scale_color_manual(values = c("top" = "#984EA3", "bottom" = "#E41A1C"), 
                         guide = F) +
      labs(x = "age (years)", y = "Adiposity trait value", 
           title = s)
    return (p)
  })
  names(res) <- STRATA
  pdf(paste0("plots/adjusted_slopes/tail_trajectories_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(as_tail_trajectories) <- PHENOTYPES

# Plot mean (binned) trajectories by adj slope quartile ----

as_quartile_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    rs <- adj_slopes[[p]][[s]]
    rs$q <- cut(rs$residual, quantile(rs$residual), include.lowest = T,
                labels = paste0("q", 1:4))
    df <- long_adiposity[[p]]
    df <- subset(df, df$eid %in% rs$eid)
    df$q <- rs$q[match(df$eid, rs$eid)]
    
    # Calculate mean and SE in each 5-year interval within each quartile
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    plot_df <- df %>% group_by(q, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    # Plot 
    p <- ggplot(plot_df, aes(x = age_bin, y = mean_value,
                             group = q, color = q, fill = q)) +
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - se_value, 
                      ymax = mean_value + se_value),
                  alpha = 0.2) +
      scale_fill_brewer(palette = "Set1") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
           title = s)
    return (p)
  })
  names(res) <- STRATA
  pdf(paste0("plots/adjusted_slopes/mean_trajectories_quartile_", p, ".pdf"),
      onefile = T)
  print(res)
  dev.off()
  return (res)
})
names(as_quartile_trajectories) <- PHENOTYPES

# Plot mean (binned) trajectories by gainer status ----

gainer_status_trajectories <- lapply(PHENOTYPES, function (p) {
  res <- lapply(STRATA, function (s) {
    rs <- adj_slopes[[p]][[s]]
    df <- long_adiposity[[p]]
    df <- subset(df, df$eid %in% rs$eid)
    df$gainer <- rs$gainer[match(df$eid, rs$eid)]
    
    # Calculate mean and SE in each 5-year interval within gainers/non-gainers
    age_bin_cuts <- seq(20, 80, by = 5)
    df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
    plot_df <- df %>% group_by(gainer, age_bin) %>% 
      summarise(count = n(),
                mean_value = mean(value),
                se_value = sd(value)/sqrt(count))
    
    plot_df$ancestry <- strsplit(s, "_")[[1]][[1]]
    plot_df$sexplot <- strsplit(s, "_")[[1]][[2]]
    return (plot_df)
  })
  res <- bind_rows(res)
  res <- split(res, res$ancestry)
  res_plots <- lapply(res, function (a) {
    # Plot 
    res_plots <- ggplot(a, aes(x = age_bin, y = mean_value,
                             group = gainer, color = gainer, fill = gainer)) +
      facet_wrap(~sexplot, nrow = 3) +
      geom_point() +
      geom_path() +
      geom_ribbon(aes(ymin = mean_value - se_value, 
                      ymax = mean_value + se_value),
                  alpha = 0.2) +
      scale_color_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
                         guide = F) +
      scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C"), 
                         guide = F) +
      labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
           title = unique(a$ancestry))
    return (res_plots)
  })
  pdf(paste0("plots/adjusted_slopes/mean_trajectories_gainer_status_", p, ".pdf"),
      onefile = T)
  print(res_plots)
  dev.off()
  return (res_plots)
})
names(gainer_status_trajectories) <- PHENOTYPES

