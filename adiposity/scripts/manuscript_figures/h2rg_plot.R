# Author: Samvida S. Venkatesh
# Date: 08/09/2022

library(tidyverse)
theme_set(theme_bw())

# colour palette: rose, teal, grey
custom_three_diverge <- c("#D35C79","#009593", "#666666")
names(custom_three_diverge) <- c("F", "M", "sex_comb")

# Read data ----

dat <- read.table("h2rg_intercept_slope.txt", sep = "\t", header = T,
                  stringsAsFactors = F)

dat <- dat %>% 
  mutate(sex = factor(sex, levels = c("F", "M", "sex_comb")),
         term = factor(term, levels = c("intercept", "slope", "correlation")),
         strata = factor(strata, levels = c("Weight_sex_comb", "Weight_M",
                                            "Weight_F",
                                            "BMI_sex_comb", "BMI_M", "BMI_F")),
         uci = beta + 1.96*se,
         lci = beta - 1.96*se)

h2rg_plot <- ggplot(dat, aes(x = beta, y = strata,
                group = term)) +
  geom_pointrange(aes(xmin = lci, xmax = uci,
                      shape = term, color = sex),
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = custom_three_diverge, guide = "none") +
  scale_x_continuous(limits = c(0, 1),
                     guide = guide_axis(check.overlap = TRUE)) +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank())

tiff("C:/Users/samvida/Documents/Lindgren Group/Adiposity_Primary_Care/Reports/Manuscript/figures/h2rg.tiff",
     height = 13, width = 5, units = "cm",
     res = 300)
print(h2rg_plot)
dev.off()
