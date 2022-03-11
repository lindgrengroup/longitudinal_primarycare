# Author: Samvida S. Venkatesh
# Date: 11/03/22

library(tidyverse)
library(ggpubr)
theme_set(theme_bw())

# Read files ----

eid_tte_matrix <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/eid_time_to_event_matrix.txt",
                             sep = "\t", header = T, stringsAsFactors = F)
colnames(eid_tte_matrix) <- gsub("^X", "", colnames(eid_tte_matrix))

# Disease dictionary
dictionary <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/UKB_codelists/chronological-map-phenotypes/annot_dictionary.txt",
                         sep = "\t", header = T, stringsAsFactors = F, 
                         quote = "", comment.char = "$")
DISEASES <- dictionary$phenotype[match(colnames(eid_tte_matrix)[c(-1:-3)],
                                       dictionary$unique_code)]
colnames(eid_tte_matrix)[c(-1:-3)] <- DISEASES

# Covariates (to add sex to file)
covars <- read.table("/well/lindgren-ukbb/projects/ukbb-11867/samvida/general_resources/220131_QCd_demographic_covariates.txt",
                     sep = "\t", header = T, comment.char = "$",
                     stringsAsFactors = F)

eid_tte_matrix$sex <- covars$sex[match(eid_tte_matrix$eid, covars$eid)]
eid_tte_matrix <- eid_tte_matrix %>% filter(!is.na(sex))

# Function to create density plot for age distributions given phenotype ----

sex_col_palette <- c("#F8766D", "#00BFC4")
names(sex_col_palette) <- c("F", "M")

plot_age_density <- function (col_name) {
  # Get data to plot
  to_plot <- eid_tte_matrix %>% 
    select(all_of(c("eid", "sex", col_name))) %>%
    filter(!is.na(!!as.symbol(col_name)))
  # Add plot sub-title on # of men and # of women with diagnosis
  nfem <- sum(to_plot$sex == "F")
  nmale <-  sum(to_plot$sex == "M")
  
  res <- ggplot(to_plot, aes(x = !!as.symbol(col_name),
                             colour = sex, fill = sex)) +
    geom_density(alpha = 0.2) +
    scale_colour_manual(values = sex_col_palette, guide = F) +
    scale_fill_manual(values = sex_col_palette, guide = F) +
    scale_x_continuous(limits = c(0, 100)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6)) +
    labs(x = "Age at first diagnosis (years)", 
         title = paste0(col_name, 
                        ", F = ", nfem,
                        ", M = ", nmale))
  return (res)
}

# Loop through all diseases and plot over multiple pages ----

plot_list <- lapply(c("age_at_first_record", "age_at_last_record",
                      DISEASES), function (pname) {
                        plot_age_density(pname)
                      })

to_print <- ggarrange(plotlist = plot_list, nrow = 3, ncol = 2)
ggexport(to_print, 
         filename = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/full_primary_care/plots/age_at_diagnosis_distributions.pdf")