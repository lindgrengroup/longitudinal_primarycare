# Author: Samvida S. Venkatesh
# Date: 30/11/2021

library(tidyverse)
library(paletteer)
theme_set(theme_bw())

# Read data ----

HORMONES <- read.table("hormone_list.txt")$V1
SEX_STRATA <- c("F", "M", "sex_comb")

enrichment_res <- lapply(HORMONES, function (hr) {
  all_hr <- lapply(SEX_STRATA, function (sx) {
    df <- read.table(paste0("results/", hr, "/phenome_wide_enrichment_", 
                            hr, "_", sx, ".txt"),
                     sep = "\t", header = T)
    df$sex_strata <- sx
    df$hormone <- hr
    return (df)
  })
  all_hr <- bind_rows(all_hr)
  return (all_hr)
})
enrichment_res <- bind_rows(enrichment_res)

annot_dictionary <- 
  read.table("C:/Users/samvida/Documents/Lindgren Group/Resources/UKBIOBANK/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt", 
             sep = "\t", header = T, stringsAsFactors = F, quote = "")

# Wrangle data for plot ----

# Highlight points that are phenome-wide significant, i.e.
# pval < 0.05 / (307 * 9 * 2) 
# This is a very conservative threshold as diseases are related to each other

# Also create a plot_y variable that is -log10(pval), but POSITIVE for over-rep
# and NEGATIVE for under-rep in order to show the difference; draw a line
# through 0 to show the split

# If pvalue is 0 from enrichment results, change it to lowest possible pval
LOWEST_P <- min(enrichment_res$pval[enrichment_res$pval != 0])
PTHRESH <- 0.05 / (307*9*2)

plot_res <- enrichment_res %>% 
  mutate(bottom_pval = ifelse(pval == 0, LOWEST_P, pval),
         highlight = bottom_pval < PTHRESH, 
         plot_y = ifelse(effect_case == "over-represented",
                         -log10(bottom_pval), log10(bottom_pval)),
         disease = factor(disease),
         hormone = factor(hormone),
         ICD_chapter = 
           annot_dictionary$ICD_chapter[match(disease, 
                                              annot_dictionary$phenotype)])

ALL_CHAPS <- unique(annot_dictionary$ICD_chapter)

col_values <- 
  colorRampPalette(paletteer_d("jcolors::pal8"))(length(ALL_CHAPS))
names(col_values) <- ALL_CHAPS

# PheWAS plot ----

all_res_phewas <- ggplot(data = plot_res %>% filter(!highlight),
                         aes(x = hormone, y = plot_y, shape = sex_strata)) +
  geom_jitter(colour = "grey", size = 1) +
  geom_jitter(data = plot_res %>% filter(highlight),
             aes(color = ICD_chapter), size = 2) +
  scale_color_manual(values = col_values, guide = "none") + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log10(PTHRESH), linetype = "dashed") +
  geom_hline(yintercept = log10(PTHRESH), linetype = "dashed") + 
  scale_y_continuous(limits = c(-200, 350), 
                     breaks = seq(-200, 350, by = 50),
                     labels = abs(seq(-200, 350, by = 50))) +
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        legend.position = "none")

ggsave("plots/phenotype_ascertainment/all_hormones.png",
    all_res_phewas,
    device = "png",
    dpi = 300, width = 15, height = 10, units = "cm")
