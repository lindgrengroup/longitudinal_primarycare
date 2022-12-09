# Author: Samvida S. Venkatesh
# Date: 30/11/2021

library(tidyverse)
library(paletteer)
theme_set(theme_bw())

# Read ad wrangle data ----

PHENOTYPES <- "Weight"
SEX_STRATA <- c("F", "M", "sex_comb")

enrichment_res <- lapply(PHENOTYPES, function (p) {
  all_pheno <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("gp_only/results/phenotype_enrichment/cluster_disease_enrichment_", 
                            p, "_", sx, ".txt"),
                     sep = "\t", header = T, quote = "", comment.char = "@")
    df <- df %>% pivot_longer(cols = starts_with("cluster"), 
                              names_to = c("cluster", ".value"),
                              names_pattern = "cluster(.)_(.*)") %>%
      mutate(neglogp = -log10(pval))
    df$sex_strata <- sx
    df$phenotype <- p
    df$disease <- gsub('"', "", df$disease)
    return (df)
  })
  names(all_pheno) <- SEX_STRATA
  return (all_pheno)
})
names(enrichment_res) <- PHENOTYPES

annot_dictionary <- 
  read.table("C:/Users/samvida/Documents/Lindgren Group/Resources/UKBIOBANK/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt", 
             sep = "\t", header = T, stringsAsFactors = F, quote = "")
annot_dictionary$phenotype <- gsub('"', "", annot_dictionary$phenotype)

# Wrangle data for plot ----

# Highlight points that are phenome-wide significant, i.e.
# pval < 0.05 / (302 * 4 * 3) 
# This is a very conservative threshold as diseases are related to each other

# Also create a plot_y variable that is -log10(pval), but POSITIVE for over-rep
# and NEGATIVE for under-rep in order to show the difference; draw a line
# through 0 to show the split

plot_res <- lapply(enrichment_res, function (sx_list) {
  res <- lapply(sx_list, function (df) {
    # If pvalue is 0 from enrichment results, change it to lowest possible pval
    LOWEST_P <- min(df$pval[df$pval != 0])
    PTHRESH <- 0.05 / (302 * 4 * 3) 
    
    for_plot <- df %>% 
      mutate(bottom_pval = ifelse(pval == 0, LOWEST_P, pval),
             highlight = bottom_pval < PTHRESH, 
             plot_y = ifelse(effect_case == "over-represented",
                             -log10(bottom_pval), log10(bottom_pval)),
             disease = factor(disease),
             cluster = factor(cluster),
             ICD_chapter = 
               annot_dictionary$ICD_chapter[match(disease, 
                                                  annot_dictionary$phenotype)])
    return (for_plot)
  })
  return (res)
})

# Create color palette
ALL_CHAPS <- sort(unique(annot_dictionary$ICD_chapter))
col_values <- 
  colorRampPalette(paletteer_d("jcolors::pal8"))(length(ALL_CHAPS))
names(col_values) <- ALL_CHAPS

# PheWAS plot for all clusters in each strata ----

all_res_phewas <- lapply(PHENOTYPES, function (p) {
  plots <- lapply(SEX_STRATA, function (sx) {
    plot_dat <- plot_res[[p]][[sx]]
    res <- ggplot(data = plot_dat %>% filter(!highlight),
                  aes(x = cluster, y = plot_y)) +
      geom_jitter(colour = "grey", size = 2) +
      geom_jitter(data = plot_dat %>% filter(highlight),
                  aes(color = ICD_chapter), size = 2) +
      scale_color_manual(values = col_values, guide = "none") + 
      geom_hline(yintercept = 0) +
      geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
      geom_hline(yintercept = log10(PTHRESH), linetype = "dashed") + 
      scale_y_continuous(limits = c(-350, 350), 
                         breaks = seq(-350, 350, by = 50),
                         labels = abs(seq(-350, 350, by = 50))) +
      theme(axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 8),
            legend.position = "none")
    
    ggsave(paste0("gp_only/plots/phenotype_enrichment/", p, "/all_clusters_", 
                  p, "_", sx, ".png"),
           res,
           device = "png",
           dpi = 300, width = 15, height = 10, units = "cm")
    return (res)
  })
  names(plots) <- SEX_STRATA
})
names(all_res_phewas) <- PHENOTYPES

# Cluster-specific plots ----

dis_levels <- annot_dictionary$phenotype[order(annot_dictionary$ICD_chapter)]

cluster_phewas <- lapply(PHENOTYPES, function (p) {
  plots <- lapply(SEX_STRATA, function (sx) {
    plot_dat <- plot_res[[p]][[sx]]
    plot_dat$disease <- factor(as.character(plot_dat$disease),
                               levels = dis_levels)
    NCLUST <- length(unique(plot_dat$cluster))
    per_clust <- lapply(1:NCLUST, function (k) {
      sub_res <- plot_dat %>% filter(cluster == k)
      lower_lim_y <- plyr::round_any(min(sub_res$plot_y), 10, f = floor)
      upper_lim_y <- plyr::round_any(max(sub_res$plot_y), 10, f = ceiling)
      breaks_y <- plyr::round_any(seq(lower_lim_y, upper_lim_y, length.out = 10),
                                  5)
      
      res <- ggplot(data = sub_res %>% filter(!highlight),
                         aes(x = disease, y = plot_y)) +
        geom_point(colour = "grey", size = 2) +
        geom_point(data = sub_res %>% filter(highlight),
                   aes(color = ICD_chapter), size = 2) +
        scale_color_manual(values = col_values, guide = "none") + 
        geom_hline(yintercept = 0) +
        geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
        geom_hline(yintercept = log10(PTHRESH), linetype = "dashed") + 
        scale_y_continuous(limits = c(lower_lim_y, upper_lim_y), 
                           breaks = breaks_y,
                           labels = abs(breaks_y)) +
        scale_x_discrete(limits = dis_levels) + 
        theme(axis.title = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 8),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              legend.position = "none")
      
      ggsave(paste0("gp_only/plots/phenotype_enrichment/", p, 
                    "/", p, "_", sx, "_cluster_", k, ".png"),
             res,
             device = "png",
             dpi = 300, width = 15, height = 10, units = "cm")
      return (res)
    })
  })
})
