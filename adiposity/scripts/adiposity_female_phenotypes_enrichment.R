# Author: Samvida S. Venkatesh
# Date: 14/04/21

library(tidyverse)
theme_set(theme_bw())
library(RColorBrewer)

# Read data ----

# Final slopes
final_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")
PHENOTYPES <- names(final_slopes)
STRATA <- names(final_slopes[[1]])
WHITE_STRATA <- STRATA[grep("white", STRATA)]

KEEP_EIDS <- lapply(PHENOTYPES, function (p) {
  eids <- lapply(WHITE_STRATA, function (s) {
    unique(final_slopes[[p]][[s]]$eid)
  })
  return (unique(unlist(eids)))
})
KEEP_EIDS <- unique(unlist(KEEP_EIDS))

# UKB phenotypes
pheno <- read.table("/well/lindgren/UKBIOBANK/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844.csv",
                    header = T, sep = ",", na.string = c("NA", "", "."), 
                    stringsAsFactors = F)
colnames(pheno) <- gsub("X", "f.", colnames(pheno))
colnames(pheno)[1] <- "eid"

pheno <- subset(pheno, pheno$eid %in% KEEP_EIDS)

# Subset columns of interest ----

# 1a. Past smoker: f.1249.0.0, f.1249.1.0, f.1249.2.0 (0 - No, 1 & 2 - Yes)
# 1b. Current smoker: f.1239.0.0, f.1239.1.0, f.1239.2.0 (1 & 2 - Yes, 3 & 4 - No)
# 2. Number of 
#   a. children fathered - f.2405.0.0, f.2405.1.0, f.2405.2.0
#   b. live births - f.2734.0.0, f.2734.1.0, f.2734.2.0
# 3. Ever used hormone replacement therapy - f.2814.0.0, f.2814.1.0, f.2814.2.0 
# (0 - No, 1 - Yes)
# 4. Recalled age at:
#   a. menarche - f.2714.0.0, f.2714.1.0, f.2714.2.0 
#   b. menopause - f.3581.0.0, f.3581.1.0, f.3581.2.0
#   c. starting HRT - f.3536.0.0, f.3536.1.0, f.3536.2.0
#   d. ending HRT - f.3546.0.0, f.3546.1.0, f.3546.2.0
# 5. Ever had hysterectomy: f.3581.0.0, f.3581.1.0, f.3581.2.0

cleaned <- pheno[, c("eid",
                     "f.1249.0.0", "f.1249.1.0", "f.1249.2.0",
                     "f.1239.0.0", "f.1239.1.0", "f.1239.2.0",
                     "f.2405.0.0", "f.2405.1.0", "f.2405.2.0",
                     "f.2734.0.0", "f.2734.1.0", "f.2734.2.0",
                     "f.2814.0.0", "f.2814.1.0", "f.2814.2.0",
                     "f.2714.0.0", "f.2714.1.0", "f.2714.2.0",
                     "f.3581.0.0", "f.3581.1.0", "f.3581.2.0",
                     "f.3536.0.0", "f.3536.1.0", "f.3536.2.0",
                     "f.3546.0.0", "f.3546.1.0", "f.3546.2.0")]
colnames(cleaned) <- c("eid", 
                       "past_smoker_1", "past_smoker_2", "past_smoker_3",
                       "current_smoker_1", "current_smoker_2", "current_smoker_3",
                       "children_fathered_1", "children_fathered_2", "children_fathered_3", 
                       "live_births_1", "live_births_2", "live_births_3", 
                       "HRT_1", "HRT_2", "HRT_3",
                       "age_menarche_1", "age_menarche_2", "age_menarche_3", 
                       "age_menopause_1", "age_menopause_2", "age_menopause_3",
                       "age_started_HRT_1", "age_started_HRT_2", "age_started_HRT_3",
                       "age_stopped_HRT_1", "age_stopped_HRT_2", "age_stopped_HRT_3")

# Create columns of interest by combining the columns above

dat <- data.frame(eid = cleaned$eid)

dat$past_smoker <- apply(cleaned %>% select(starts_with("past_smoker")),
                         1, function (r) {
                           ifelse(any(r %in% c(1, 2), na.rm = T), T,
                                  ifelse(any(r %in% c(3, 4), na.rm = T), F, 
                                         NA))
                         })

dat$current_smoker <- apply(cleaned %>% select(starts_with("current_smoker")),
                            1, function (r) {
                              ifelse(any(r %in% c(1, 2), na.rm = T), T,
                                     ifelse(any(r == 0, na.rm = T), F, 
                                            NA))
                            })

dat$n_children_fathered <- apply(cleaned %>% select(starts_with("children")),
                                 1, function (r) { 
                                   ifelse(all(is.na(r)), NA, 
                                          ifelse(any(r < 0, na.rm = T), NA,
                                                 median(r, na.rm = T)))
                                 })

dat$n_live_births <- apply(cleaned %>% select(starts_with("live")),
                           1, function (r) { 
                             ifelse(all(is.na(r)), NA, 
                                    ifelse(any(r < 0, na.rm = T), NA,
                                           median(r, na.rm = T)))
                           })

dat$ever_HRT <- apply(cleaned %>% select(starts_with("HRT")),
                      1, function (r) {
                        ifelse(any(r == 1, na.rm = T), T,
                               ifelse(any(r == 0, na.rm = T), F, 
                                      NA))
                      })

dat$age_menarche <- apply(cleaned %>% select(starts_with("age_menarche")),
                          1, function (r) { 
                            ifelse(all(is.na(r)), NA, 
                                   ifelse(any(r < 0, na.rm = T), NA,
                                          median(r, na.rm = T)))
                          })

dat$age_menopause <- apply(cleaned %>% select(starts_with("age_menopause")),
                           1, function (r) { 
                             ifelse(all(is.na(r)), NA, 
                                    ifelse(any(r < 0, na.rm = T), NA,
                                           median(r, na.rm = T)))
                           })

dat$age_started_HRT <- apply(cleaned %>% select(starts_with("age_started_HRT")),
                             1, function (r) { 
                               ifelse(all(is.na(r)), NA, 
                                      ifelse(any(r < 0, na.rm = T), NA,
                                             median(r, na.rm = T)))
                             })

dat$age_stopped_HRT <- apply(cleaned %>% select(starts_with("age_stopped_HRT")),
                             1, function (r) { 
                               ifelse(all(is.na(r)), NA, 
                                      ifelse(any(r < 0, na.rm = T), NA,
                                             median(r, na.rm = T)))
                             })

# Add UKB-phenotypes to slopes dataframes ----

for_testing <- lapply(PHENOTYPES, function (p) {
  df <- lapply(WHITE_STRATA, function (s) {
    res <- merge(final_slopes[[p]][[s]][, c("eid", "sex", "residual", "gainer")],
                 dat, by = "eid")
    return (res)
  })
  names(df) <- WHITE_STRATA
  return (df)
})
names(for_testing) <- PHENOTYPES

saveRDS(for_testing, "results/phenotype_analysis/female_phenotypes.rds")

# Plot distributions of UKB-phenotypes by gainer status ----

pheno_plots <- lapply(PHENOTYPES, function (p) {
  plot_lists <- lapply(WHITE_STRATA, function (s) {
    df <- for_testing[[p]][[s]]
    
    # Ages at menarche, menopause, started HRT, ended HRT
    df_ages <- pivot_longer(df, cols = c(starts_with("age_")),
                            names_to = "event", values_to = "age")
    if (!all(is.na(df_ages$age))) {
      p1 <- ggplot(df_ages, aes(x = gainer, y = age)) +
        facet_wrap(~event, nrow = 2) +
        geom_violin(aes(fill = gainer), position = position_dodge(1)) +
        geom_boxplot(width = 0.1) + 
        scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C")) + 
        labs(x = "gainer status", y = "age at event", title = s)
    } else { p1 <- NA }
    
    # Proportion to report past smoking, current smoking, having had HRT
    df_props <- df[, c("eid", "gainer", "past_smoker", "current_smoker",
                       "ever_HRT")] %>% group_by(gainer) %>%
      summarise_at(c("past_smoker", "current_smoker", "ever_HRT"), 
                   mean, na.rm = T)
    df_props <- pivot_longer(df_props, cols = -gainer,
                             names_to = "attribute", values_to = "percentage")
    
    p2 <- ggplot(df_props, aes(x = attribute, y = percentage,
                               color = gainer, fill = gainer)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C")) + 
      scale_color_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C")) + 
      #scale_x_discrete(labels = c("past smoker", "current smoker",
      #"ever used HRT")) +
      labs(y = "Frequency in cohort", title = s)
    
    # Number of live births and n children fathered
    df_n <- pivot_longer(df, cols = c(starts_with("n_")),
                         names_to = "event", values_to = "number") %>%
      mutate(number = floor(number))
    
    p3 <- ggplot(df_n, aes(x = gainer, y = number)) +
      facet_wrap(~event, ncol = 2) +
      geom_violin(aes(fill = gainer), position = position_dodge(1)) +
      geom_boxplot(width = 0.1) + 
      scale_fill_manual(values = c("TRUE" = "#984EA3", "FALSE" = "#E41A1C")) + 
      labs(x = "gainer status", y = "number of children (fathered/birthed)",
           title = s)
    
    return (list(p1 = p1, p2 = p2, p3 = p3))
    
  })
  
  names(plot_lists) <- WHITE_STRATA
  
  pdf(paste0("plots/phenotype_analysis/female_phenotypes_compare_", p, ".pdf"), 
      onefile = T)
  lapply(plot_lists, function (res) {
    print(res)
  })
  dev.off()
  
  return(plot_lists)
})
names(pheno_plots) <- PHENOTYPES

# Statistical tests for differences ----

test_results <- lapply(PHENOTYPES, function (p) {
  res <- lapply(WHITE_STRATA, function (s) {
    df <- for_testing[[p]][[s]]
    
    # Rename gainer vs non-gainer columns
    df$gainer <- ifelse(df$gainer, "gainer", "non_gainer")
    # t-test for continuous variables:
    # Ages at menarche, menopause, started HRT, ended HRT,
    # n-chidren-fathered, n-live-births
    if (s == "white_M") {
      t_test_res <- df %>% 
        select(gainer, n_children_fathered) %>% 
        gather(key = variable, value = value, -gainer) %>% 
        group_by(gainer, variable) %>% 
        summarise(value = list(value)) %>% 
        spread(gainer, value) %>% 
        group_by(variable) %>% 
        mutate(p_value = t.test(unlist(gainer), unlist(non_gainer))$p.value,
               t_value = t.test(unlist(gainer), unlist(non_gainer))$statistic)
    } else if (s == "white_F") {
      t_test_res <- df %>% 
        select(gainer, age_menarche, age_menopause, 
               age_started_HRT, age_stopped_HRT, 
               n_live_births) %>% 
        gather(key = variable, value = value, -gainer) %>% 
        group_by(gainer, variable) %>% 
        summarise(value = list(value)) %>% 
        spread(gainer, value) %>% 
        group_by(variable) %>% 
        mutate(p_value = t.test(unlist(gainer), unlist(non_gainer))$p.value,
               t_value = t.test(unlist(gainer), unlist(non_gainer))$statistic)
    } else {
      t_test_res <- df %>% 
        select(gainer, age_menarche, age_menopause, 
               age_started_HRT, age_stopped_HRT, 
               n_live_births, n_children_fathered) %>% 
        gather(key = variable, value = value, -gainer) %>% 
        group_by(gainer, variable) %>% 
        summarise(value = list(value)) %>% 
        spread(gainer, value) %>% 
        group_by(variable) %>% 
        mutate(p_value = t.test(unlist(gainer), unlist(non_gainer))$p.value,
               t_value = t.test(unlist(gainer), unlist(non_gainer))$statistic)
    }
    t_test_res <- t_test_res[, c("variable", "p_value", "t_value")]
    
    write.table(t_test_res, 
                paste0("results/phenotype_analysis/female_phenotypes_continuous_",
                p, "_", s, ".txt"), sep = "\t", row.names = F, quote = F)
    
    # Fisher's test for binary variables
    # past smoker, current smoker, and ever_HRT
    
    # Number to report past smoking, having had HRT
    df_props <- df[, c("eid", "gainer", "past_smoker", "ever_HRT")] %>% 
      group_by(gainer) %>%
      summarise_at(c("past_smoker", "ever_HRT"), sum, na.rm = T)
    NGAINERS <- sum(df$gainer == "gainer")
    NCOHORT <- nrow(df)
    # gainers vs non-gainers
    # past smoking
    mat_past_smoker <- matrix(c(df_props$past_smoker[df_props$gainer == "gainer"],
                      sum(df_props$past_smoker),
                       NGAINERS - df_props$past_smoker[df_props$gainer == "gainer"],
                       NCOHORT - sum(df_props$past_smoker)), 2, 2)
    pval_past_smoker <- fisher.test(mat_past_smoker)$p.value
    
    if (s == "white_M") {
      pval_HRT <- NA
    } else {
      mat_HRT <- matrix(c(df_props$ever_HRT[df_props$gainer == "gainer"],
                                  sum(df_props$ever_HRT),
                                  NGAINERS - df_props$ever_HRT[df_props$gainer == "gainer"],
                                  NCOHORT - sum(df_props$ever_HRT)), 2, 2)
      pval_HRT <- fisher.test(mat_HRT)$p.value
    }
    
    sink(paste0("results/phenotype_analysis/female_phenotypes_binary_",
                p, "_", s, ".txt"))
    cat("pval_past_smoker: ", pval_past_smoker, "\n")
    cat("pval_HRT: ", pval_HRT, "\n")
    sink()
    
    return (list(t_test = t_test_res, 
                 past_smoking = pval_past_smoker, 
                 HRT = pval_HRT))
  })
  
  names(res) <- WHITE_STRATA
  return(res)
})
names(test_results) <- PHENOTYPES
