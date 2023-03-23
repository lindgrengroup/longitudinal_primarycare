# Author: Samvida S. Venkatesh
# Date: 11/10/21

library(tidyverse)
library(ggrepel)
library(paletteer)
theme_set(theme_bw())

# Read data ----

mainpath <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/binary_phewas" # REDACTED
gen_resources_path <- "/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources" # REDACTED

res <- read.table(paste0(mainpath, "/rs429358_all_results.txt"), 
                  sep = "\t",
                  header = T, stringsAsFactors = F, quote = "", comment.char = "@")
# Disease dictionary
dictionary <- read.table(paste0(gen_resources_path, "/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt"),
                         sep = "\t", header = T, stringsAsFactors = F,
                         comment.char = "@", quote = "")
# colour palette
ALL_CHAPS <- sort(unique(dictionary$ICD_chapter))

col_values <- 
  colorRampPalette(paletteer_d("jcolors::pal8"))(length(ALL_CHAPS))
names(col_values) <- ALL_CHAPS

# Wrangle data for plot ----

# Highlight points that are phenome-wide significant, i.e.
# pval < 0.05 / DX_LEVELS
# This is a very conservative threshold as diseases are related to each other

# Also create a plot_y variable that is -log10(pval), but POSITIVE for over-rep
# and NEGATIVE for under-rep in order to show the difference; draw a line
# through 0 to show the split

alpha_vals <- c(0.5, 1)
names(alpha_vals) <- c("no", "signif")

# If pvalue is 0 from enrichment results, change it to lowest possible pval
LOWEST_P <- min(res$pvalue[res$pvalue != 0])
PTHRESH <- 0.05 / length(unique(res$disease))

plot_dat <- res %>% 
  mutate(bottom_pval = ifelse(pvalue == 0, LOWEST_P, pvalue),
         highlight = ifelse(bottom_pval <= PTHRESH, "signif", "no"),
         plot_y = ifelse(OR >= 1,
                         -log10(bottom_pval), log10(bottom_pval)),
         ICD_chapter = 
           dictionary$ICD_chapter[match(disease, 
                                        dictionary$phenotype)],
         ICD_chapter = factor(ICD_chapter, levels = ALL_CHAPS)) %>%
  group_by(ICD_chapter) %>%
  arrange(disease, .by_group = T)

DX_ORDER <- unique(plot_dat$disease)
plot_dat$disease <- factor(plot_dat$disease, levels = DX_ORDER)

getPlot <- function (df) {
  lower_lim_y <- plyr::round_any(min(df$plot_y), 10, f = floor)
  upper_lim_y <- plyr::round_any(max(df$plot_y), 10, f = ceiling)
  breaks_y <- plyr::round_any(seq(lower_lim_y, upper_lim_y, length.out = 10),
                              5)
  
  res_plot <- ggplot(data = df,
                     aes(x = disease, y = plot_y, 
                         color = ICD_chapter, fill = ICD_chapter)) +
    geom_point(aes(alpha = highlight)) +
    geom_text_repel(data = df %>% filter(highlight == "signif"), 
                    aes(label = disease), size = 3,
                    show.legend = F, box.padding = 0.5, max.overlaps = Inf) + 
    scale_alpha_manual(values = alpha_vals, guide = "none") + 
    scale_fill_manual(values = col_values) +
    scale_colour_manual(values = col_values) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
    geom_hline(yintercept = log10(PTHRESH), linetype = "dashed") + 
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

SEX_STRATA <- c("F", "M", "sex_comb")
lapply(SEX_STRATA, function (sx) {
  sub_res <- plot_dat %>% filter(sex_strata == sx &
                                   abs(plot_y) < 100)
  
  ggsave(paste0(mainpath, "/plots/phewas_", sx, ".png"),
         getPlot(sub_res),
         width = 6.5, height = 2.5, units = "in")
})
