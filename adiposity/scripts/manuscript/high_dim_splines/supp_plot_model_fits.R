# Author: Samvida S. Venkatesh
# Adapted from: George Nicholson
# Date: 16/05/22

library(argparse)
library(splines)
library(tidyverse)
library(RColorBrewer)
theme_set(theme_bw())

RANDOM_SEED <- 160522
set.seed(RANDOM_SEED)

mainpath <- "" # REDACTED
gen_resources_path <- "" # REDACTED
hidim_mods_path <- "" # REDACTED

parser <- ArgumentParser()
parser$add_argument("--phenotype", required = TRUE,
                    help = "Phenotype to model")
parser$add_argument("--sex_strata", required = TRUE,
                    help = "Sex strata")
args <- parser$parse_args()

PHENO <- args$phenotype
SEX_STRATA <- args$sex_strata

plotdir <- paste0(hidim_mods_path, "/plots/")
resdir <- paste0(hidim_mods_path, "/results/")

# Load data ----

raw_dat <- readRDS(paste0(hidim_mods_path, "/data/dat_to_model.rds"))[[PHENO]][[SEX_STRATA]]
model_dat <- readRDS(paste0(hidim_mods_path, "/results/fit_objects_", 
                            PHENO, "_", SEX_STRATA, ".rds"))

B <- model_dat$B
spline_posteriors <- model_dat$spline_posteriors
model_resid_var <- model_dat$resid_var

# Plot fitted values for sample data ----

## only run this after getting the residual variance from above

set.seed(RANDOM_SEED)
plot_indivs <- sample(unique(raw_dat$eid), 25, replace = F)

pred_sample <- lapply(plot_indivs, function (id) {
  pred_value <- B %*% spline_posteriors[[id]]$mu
  sd_pred <- sqrt(diag(B %*% spline_posteriors[[id]]$Sig %*% t(B)) * model_resid_var)

  res <- data.frame(eid = id,
                    t_diff = 1:length(pred_value),
                    fit_mean = pred_value,
                    fit_sd = sd_pred)
  return (res)
})
pred_sample <- bind_rows(pred_sample)

raw_sample <- raw_dat %>% filter(eid %in% plot_indivs)

fit_plot <- ggplot(pred_sample, aes(x = t_diff)) +
  facet_wrap(~eid, nrow = 5, ncol = 5, scales = "free_y") +
  geom_point(data = raw_sample,
             aes(x = t_diff, y = value_fulladj)) +
  geom_line(aes(y = fit_mean)) +
  geom_ribbon(aes(ymin = fit_mean - 1.96*fit_sd,
                  ymax = fit_mean + 1.96*fit_sd), alpha = 0.1) +
  labs(x = "Days from first measurement", y = "Confounder-adj value")

pdf(paste0(plotdir, "sample_pred_", PHENO, "_", SEX_STRATA, ".pdf"), onefile = T)
print(fit_plot)
dev.off()
