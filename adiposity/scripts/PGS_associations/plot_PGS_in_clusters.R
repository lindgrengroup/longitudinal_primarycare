# Author: Samvida S. Venkatesh
# Date: 08/02/2022

library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

# Read files ----

PHENOTYPES <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/pheno_names.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")
PGS_LIST <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/ukb_pgs_directories_list.txt",
                       sep = "\t", header = T, stringsAsFactors = F)
PGS_TRAITS <- PGS_LIST$trait

general_covars <- read.table("/well/lindgren/UKBIOBANK/samvida/general_resources/QCd_demographic_covariates.txt",
                             sep = "\t", header = T, comment.char = "$",
                             stringsAsFactors = F)

trait_covars <- readRDS("/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/data/covariates.rds")[PHENOTYPES]

cluster_membership <- lapply(PHENOTYPES, function (p) {
  res_list <- lapply(SEX_STRATA, function (sx) {
    res <- read.table(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/results/cubic_splines/not_adj_genetic_PCs/time_based/age_adj_clusters_",
                      p, "_", sx, ".txt"),
               sep = "\t", header = T, stringsAsFactors = F)
    res$eid <- as.character(res$eid)
    res$baseline_trait <- 
      trait_covars[[p]]$baseline_trait[match(res$eid, trait_covars[[p]]$eid)]
    return (res)
  })
  names(res_list) <- SEX_STRATA
  return (res_list)
})
names(cluster_membership) <- PHENOTYPES

pgs <- lapply(1:nrow(PGS_LIST), function (i) {
  df <- read.table(PGS_LIST$pgs_location[i], header = F, stringsAsFactors = F)
  colnames(df) <- c("eid", paste0("PGS_", PGS_TRAITS[i]))
  df$eid <- as.character(df$eid)
  return (df)
})
pgs <- pgs %>% reduce(full_join, by = "eid")
pgs$sex <- general_covars$sex[match(pgs$eid, general_covars$eid)]

# Wrangle data for plot ----

dat <- lapply(PHENOTYPES, function (p) {
  dat_list <- lapply(SEX_STRATA, function (sx) {
    
    clust_dat <- cluster_membership[[p]][[sx]] %>%
      mutate(k = as.character(k), 
             ktype = "cubic_spline_clusters")
    
    # Add baseline trait quantiles as cluster of different type
    ncluster <- length(unique(clust_dat$k))
    quant_dat <- clust_dat %>% 
      mutate(k = as.character(ntile(baseline_trait, ncluster)),
             ktype = "baseline_trait_quantiles")  
    
    plot_dat <- bind_rows(clust_dat, quant_dat)
    plot_dat <- left_join(plot_dat, pgs, by = "eid")
    
    # Add "full UKB" as cluster and include all individuals of the
    # relevant sex
    add_rows <- pgs
    if (sx != "sex_comb") add_rows <- add_rows %>% filter(sex == sx)
    add_rows <- add_rows %>% 
      select(-any_of("sex")) %>% 
      mutate(k = "full_UKB")
    
    add_rows1 <- add_rows %>% mutate(ktype = "cubic_spline_clusters")
    add_rows2 <- add_rows %>% mutate(ktype = "baseline_trait_quantiles")  
    
    plot_dat <- bind_rows(plot_dat, add_rows1, add_rows2)
    plot_dat$k <- factor(as.character(plot_dat$k))
    return (plot_dat)
  })
  names(dat_list) <- SEX_STRATA
  return (dat_list)
})
names(dat) <- PHENOTYPES

# Plot cluster identity vs PGSs ----

plot_cluster_vs_pgs <- function (p, sx, pgs_trait) {
  pgs_colname <- paste0("PGS_", pgs_trait)
   pgs_plot <- ggplot(dat[[p]][[sx]], 
                        aes(x = k, y = get(pgs_colname))) +
    facet_wrap(~ktype, nrow = 2) +
    geom_violin(aes(fill = k), position = position_dodge(1)) +
    geom_boxplot(width = 0.1) + 
    scale_fill_brewer(palette = "Dark2") + 
    labs(x = "Cluster identity", y = pgs_colname, 
         title = paste0(pgs_colname, " distribution in clusters of: ",
                        p, " strata: ", sx)) +
    theme(legend.position = "none")

  return (pgs_plot)
}

# Apply all plotting functions ----

lapply(PHENOTYPES, function (p) {
  lapply(SEX_STRATA, function (sx) {
    per_pgs_trait <- lapply(PGS_TRAITS, function (pt) {
      plot_cluster_vs_pgs(p, sx, pt)
    })
    pdf(paste0("/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/gp_only/plots/cubic_splines/clustering/cluster_pgs_associations_",
               p, "_", sx, ".pdf"))
    print(per_pgs_trait)
    dev.off()
  })
})
