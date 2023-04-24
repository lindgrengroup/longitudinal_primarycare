# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)
library(ggrepel)
library(paletteer)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/binary_phewas" # REDACTED
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources" # REDACTED

SEX_STRATA <- c("F", "M", "sex_comb")

VARIDS <- c("6_26076446_GA_G", "rs429358","rs1175592", "rs9467663", 
            "rs12953815", "rs61955499")
# Should the ORs plotted be flipped to align with the BMI/weight-increasing allele?
flip_to_incr <- c(F, T, F, T, T, T)
names(flip_to_incr) <- VARIDS

all_res <- lapply(VARIDS, function (varid) {
  res <- read.table(paste0(mainpath, "/", varid, "_all_results.txt"), 
                    sep = "\t",
                    header = T, stringsAsFactors = F, quote = "", comment.char = "@")
  return (res)
})
names(all_res) <- VARIDS

# Disease dictionary
dictionary <- read.table(paste0(gen_resources_path, "/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt"),
                         sep = "\t", header = T, stringsAsFactors = F,
                         comment.char = "@", quote = "")
# colour palette
ALL_CHAPS <- sort(unique(dictionary$ICD_chapter))

col_values <- 
  colorRampPalette(paletteer_d("jcolors::pal8"))(length(ALL_CHAPS))
names(col_values) <- ALL_CHAPS

alpha_vals <- c(0.5, 1)
names(alpha_vals) <- c("no", "signif")

# pval < 0.05 / DX_LEVELS * 6 
PTHRESH <- 0.05 / (length(unique(dictionary$phenotype)) * 6)

# Wrangle data for plot ----

# Highlight points that are phenome-wide significant, i.e.
# This is a very conservative threshold as diseases are related to each other

# Also create a plot_y variable that is -log10(pval), but POSITIVE for over-rep
# and NEGATIVE for under-rep in order to show the difference; draw a line
# through 0 to show the split

wrangleDat <- function (raw_df) {
  # If pvalue is 0 from enrichment results, change it to lowest possible pval
  lowest_p <- min(raw_df$pvalue[raw_df$pvalue != 0])
  
  plot_dat <- raw_df %>% 
    mutate(bottom_pval = ifelse(pvalue == 0, lowest_p, pvalue),
           highlight = ifelse(bottom_pval <= PTHRESH, "signif", "no"),
           plot_y = ifelse(plot_OR >= 1,
                           -log10(bottom_pval), log10(bottom_pval)),
           ICD_chapter = 
             dictionary$ICD_chapter[match(disease, 
                                          dictionary$phenotype)],
           ICD_chapter = factor(ICD_chapter, levels = ALL_CHAPS)) %>%
    group_by(ICD_chapter) %>%
    arrange(disease, .by_group = T)
  
  dx_order <- unique(plot_dat$disease)
  plot_dat$disease <- factor(plot_dat$disease, levels = dx_order)
  return (plot_dat)
}

getPlot <- function (df) {
  lower_lim_y <- plyr::round_any(min(df$plot_y), 10, f = floor)
  upper_lim_y <- plyr::round_any(max(df$plot_y), 10, f = ceiling)
  breaks_y <- plyr::round_any(seq(lower_lim_y, upper_lim_y, length.out = 10),
                              5)
  
  res_plot <- ggplot(data = df,
                     aes(x = disease, y = plot_y, 
                         color = ICD_chapter, fill = ICD_chapter)) +
    geom_point(aes(alpha = highlight)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
    geom_hline(yintercept = log10(PTHRESH), linetype = "dashed") + 
    geom_text_repel(data = df %>% filter(highlight == "signif"), 
                    aes(label = disease), size = 3,
                    show.legend = F, box.padding = 0.5, max.overlaps = Inf) + 
    scale_alpha_manual(values = alpha_vals, guide = "none") + 
    scale_fill_manual(values = col_values) +
    scale_colour_manual(values = col_values) +
    scale_y_continuous(limits = c(lower_lim_y, upper_lim_y), 
                       breaks = breaks_y,
                       labels = abs(breaks_y)) +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  return (res_plot)
}

lapply(VARIDS, function (varid) {
  # Prepare data
  raw_df <- all_res[[varid]]
  raw_df$plot_OR <- raw_df$OR
  if (flip_to_incr[[varid]]) {
    raw_df$plot_OR <- 1/raw_df$OR
  }
  
  lapply(SEX_STRATA, function (sx) {
    wrangled_df <- wrangleDat(raw_df)
    sub_res <- wrangled_df %>% filter(sex_strata == sx &
                                     abs(plot_y) < 100)
    
    ggsave(paste0(mainpath, "/plots/", varid, "_", sx, ".png"),
           getPlot(sub_res),
           width = 6.5, height = 2.5, units = "in")
  })
})

