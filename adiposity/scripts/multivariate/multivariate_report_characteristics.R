# Author: Samvida S. Venkatesh
# Date: 19/03/21

library(tidyverse)
library(RColorBrewer)
library(UpSetR)
theme_set(theme_bw())

# Read stratified slope and adiposity data ----

adj_slopes <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/adiposity_adj_slopes.rds")
adiposity <- readRDS("/well/lindgren/UKBIOBANK/samvida/adiposity/stratified_adiposity.rds")
PHENOTYPES <- names(adj_slopes)
STRATA <- unique(unlist(lapply(adj_slopes, function (x) names(x) )))

SC_STRATA <- STRATA[grepl("sexcomb", STRATA)]

# Get list of individuals with all adiposity traits ----

adipo_ids <- lapply(PHENOTYPES, function (p) {
  id_vec <- lapply(SC_STRATA, function (s) { adj_slopes[[p]][[s]]$eid })
  return (unlist(id_vec))
})
# Get intersection of IDs that are present in all adiposity datasets
MULTI_IDS <- Reduce(intersect, adipo_ids)

# Get list of gainers by strata ----

# Create lists
strata_list <- lapply(STRATA, function (s) {
  p_gainers <- lapply(PHENOTYPES, function (p) {
    rs <- adj_slopes[[p]][[s]]
    # Retain only multi-trait individuals 
    rs <- subset(rs, rs$eid %in% MULTI_IDS)
    # Create lists of gainers for upset plot
    return (rs$eid[rs$gainers])
  })
  names(p_gainers) <- paste(PHENOTYPES, "_gainers")
  return (p_gainers)
})
names(strata_list) <- STRATA

pdf("plots/multivariate/adj_slope_gainer_comparison.pdf", onefile = T)
lapply(STRATA, function (s) {
  p <- upset(fromList(strata_list[[s]]), order.by = "freq", 
             mainbar.y.label = paste("# of", s))
  print (p)
  return()
})
dev.off()

# Plot adiposity trajectories for individuals in gainer intersections ----

# Focus on white F for these visualisations (but can extend to all)
s <- "white_F"

sub_slopes <- lapply(PHENOTYPES, function (p) {
  res <- adj_slopes[[p]][[s]]
  res <- subset(res, res$eid %in% MULTI_IDS)
  res <- res[, c("eid", "gainers")]
  colnames(res) <- c("eid", paste(p, "gainer", sep = "_"))
  return (res)
})
sub_slopes <- Reduce(merge, sub_slopes)

# Get relevant adiposity trajectories
sub_adipo <- lapply(PHENOTYPES, function (p) {
  res <- adiposity[[p]][[s]]
  res$adipo_trait <- p
  return (res)
})
sub_adipo <- bind_rows(sub_adipo)

sub_adipo$adipo_trait <- factor(sub_adipo$adipo_trait, 
                                levels = c("BMI", "weight", "waist", "WHR"))

## Get groups ----

gp1 <- list()
gp2 <- list()
sets <- c("all gainers v. none gainers",
          "BMI and weight gainers v. WC and WHR gainers", 
          "WHR gainers v. other gainers",
          "WC gainers v. other gainers",
          "weight gainers v. other gainers",
          "BMI gainers v. other gainers")

# Set 1: all gainers vs none gainers
gp1[[sets[1]]] <- sub_slopes$eid[rowSums(sub_slopes[, 2:5]) == 4]
gp2[[sets[1]]] <- sub_slopes$eid[rowSums(sub_slopes[, 2:5]) == 0]

# Set 2: BMI and weight gainers v. WC and WHR gainers
gp1[[sets[2]]] <- sub_slopes$eid[sub_slopes$BMI_gainer & 
                                   sub_slopes$weight_gainer &
                                   !sub_slopes$waist_gainer & 
                                   !sub_slopes$WHR_gainer]
gp2[[sets[2]]] <- sub_slopes$eid[!sub_slopes$BMI_gainer & 
                                   !sub_slopes$weight_gainer &
                                   sub_slopes$waist_gainer & 
                                   sub_slopes$WHR_gainer]

# Set 3: WHR gainers v. other gainers
gp1[[sets[3]]] <- sub_slopes$eid[!sub_slopes$BMI_gainer & 
                                   !sub_slopes$weight_gainer &
                                   !sub_slopes$waist_gainer & 
                                   sub_slopes$WHR_gainer]
gp2[[sets[3]]] <- sub_slopes$eid[sub_slopes$BMI_gainer & 
                                   sub_slopes$weight_gainer &
                                   sub_slopes$waist_gainer & 
                                   !sub_slopes$WHR_gainer]

# Set 4: WC gainers v. other gainers
gp1[[sets[4]]] <- sub_slopes$eid[!sub_slopes$BMI_gainer & 
                                   !sub_slopes$weight_gainer &
                                   sub_slopes$waist_gainer & 
                                   !sub_slopes$WHR_gainer]
gp2[[sets[4]]] <- sub_slopes$eid[sub_slopes$BMI_gainer & 
                                   sub_slopes$weight_gainer &
                                   !sub_slopes$waist_gainer & 
                                   sub_slopes$WHR_gainer]

# Set 5: weight gainers v. other gainers
gp1[[sets[5]]] <- sub_slopes$eid[!sub_slopes$BMI_gainer & 
                                   sub_slopes$weight_gainer &
                                   !sub_slopes$waist_gainer & 
                                   !sub_slopes$WHR_gainer]
gp2[[sets[5]]] <- sub_slopes$eid[sub_slopes$BMI_gainer & 
                                   !sub_slopes$weight_gainer &
                                   sub_slopes$waist_gainer & 
                                   sub_slopes$WHR_gainer]

# Set 6: BMI gainers v. other gainers
gp1[[sets[6]]] <- sub_slopes$eid[sub_slopes$BMI_gainer & 
                                   !sub_slopes$weight_gainer &
                                   !sub_slopes$waist_gainer & 
                                   !sub_slopes$WHR_gainer]
gp2[[sets[6]]] <- sub_slopes$eid[!sub_slopes$BMI_gainer & 
                                   sub_slopes$weight_gainer &
                                   sub_slopes$waist_gainer & 
                                   sub_slopes$WHR_gainer]

plot_dfs <- lapply(sets, function (set) {
  res <- subset(sub_adipo, 
                sub_adipo$eid %in% gp1[[set]] | sub_adipo$eid %in% gp2[[set]])
  res$category <- ifelse(res$eid %in% gp1[[set]], "gp1", "gp2")
  return (res)
})
names(plot_dfs) <- sets

## Binned mean trajectories by group ----

age_bin_cuts <- seq(20, 80, by = 5)

bin_plots <- lapply(sets, function (set) {
  df <- plot_dfs[[set]] 
  df$age_bin <- cut(df$age_event, age_bin_cuts, include.lowest = T)
  df <- df %>% group_by(adipo_trait, category, age_bin) %>% 
    summarise(count = n(),
              mean_value = mean(value),
              se_value = sd(value)/sqrt(count))
  # Plot 
  p <- ggplot(df, aes(x = age_bin, y = mean_value,
                          group = category, color = category, fill = category)) +
    facet_wrap(~adipo_trait, nrow = 2, scale = "free_y") +
    geom_point() +
    geom_path() +
    geom_ribbon(aes(ymin = mean_value - se_value, 
                    ymax = mean_value + se_value),
                alpha = 0.2) +
    scale_fill_manual(values = c("gp1" = "#984EA3", 
                                  "gp2" = "#E41A1C")) +
    scale_color_manual(values = c("gp1" = "#984EA3", 
                                  "gp2" = "#E41A1C")) +
    labs(x = "Age bin (years)", y = "Mean (S.E.) of adiposity",
         title = set)
  
  return (p)
})

pdf("plots/multivariate/white_F_multivariate_groups_binned.pdf", onefile = T)
print(bin_plots)
dev.off()

## Sample trajectories by group ----
sample_plots <- lapply(sets, function (set) {
  df <- plot_dfs[[set]] 
  # Get sample trajectories 
  plot_ids_gp1 <- sample(gp1[[set]], 1, replace = F)
  plot_ids_gp2 <- sample(gp2[[set]], 1, replace = F)
  df <- subset(df, df$eid %in% plot_ids_gp1 | df$eid %in% plot_ids_gp2)
  df$category <- ifelse(df$eid %in% plot_ids_gp1, "gp1", "gp2")
  # Plot
  p <- ggplot(df, aes(x = age_event, y = value, 
                           group = eid, color = category)) +
    facet_wrap(~adipo_trait, nrow = 2, scale = "free_y") +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c("gp1" = "#984EA3", 
                                  "gp2" = "#E41A1C")) +
    labs(x = "Age (years)", y = "Adiposity trait value", 
         title = set)
  return (p)
})

pdf("plots/multivariate/white_F_multivariate_groups_sample.pdf", onefile = T)
print(sample_plots)
dev.off()



# FOCUS: BMI AND WHR in white-ancestry individuals ----

white_strata <- STRATA[grep("white", STRATA)]
two_by_two <- lapply(white_strata, function (s) {
  WHR_df <- adj_slopes[["WHR"]][[s]]
  waist_df <- adj_slopes[["waist"]][[s]]
  
  WHRg <- WHR_df$eid[WHR_df$gainer]
  WHRl <- WHR_df$eid[!WHR_df$gainer]
  
  waistg <- waist_df$eid[waist_df$gainer]
  waistl <- waist_df$eid[!waist_df$gainer]
  
  gg = length(intersect(WHRg, waistg))
  gl = length(intersect(WHRg, waistl))
  lg = length(intersect(WHRl, waistg))
  ll = length(intersect(WHRl, waistl))
  tot <- gg + gl + lg + ll
  
  res <- data.frame(WHRg = c(gg/tot, gl/tot), WHRl = c(lg/tot, ll/tot))
  rownames(res) <- c("waistg", "waistl")
  
  return (res)
})
names(two_by_two) <- white_strata

