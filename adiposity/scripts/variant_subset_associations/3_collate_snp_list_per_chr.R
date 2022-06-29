# Author: Samvida S. Venkatesh
# Date: 29/06/22

library(tidyverse)
library(biomaRt)

main_filepath <- "/well/lindgren/samvida/Resources/"

# Read GWAS catalog variants and traits to extract ----

gwas_cat_snps <- read.table(paste0(main_filepath, "GWASCatalog/all_associations_211102.txt"),
                            sep = "\t", header = T, quote = "")
# Traits to extract

gwas_cat_obesity_traits <- read.table(paste0(main_filepath, "GWASCatalog/gwascat_obesity_traits.txt"),
                                      sep = "\t", header = F, stringsAsFactors = F)$V1

hg38_chr_pos <- read.table(paste0(main_filepath, "hg38/grch38.all.rsid.txt"),
                           sep = " ", header = T, quote = "")

# Get associations ----

obesity_snps <- gwas_cat_snps %>% 
  filter(MAPPED_TRAIT %in% gwas_cat_obesity_traits) %>%
  select(all_of(c("CHR_ID", "CHR_POS", "SNPS"))) %>%
  distinct()

# Add chr and pos if missing
obesity_snps <- obesity_snps %>% 
  mutate(CHR_add = hg38_chr_pos$chr[match(SNPS, hg38_chr_pos$snp)],
         POS_add = hg38_chr_pos$pos[match(SNPS, hg38_chr_pos$snp)]) %>%
  # flag variants with incorrect chr positions
  mutate(flag = ifelse(CHR_add != CHR_ID | POS_add != CHR_POS))

# Comment out later
write.table(obesity_snps, 
            paste0(main_filepath, "GWASCatalog/tmp_obesity_associations.txt"),
            sep = "\t", row.names = F, quote = F, col.names = T)

# Write bed file so it can be lifted over
obesity_bed <- obesity_snps[complete.cases(obesity_snps), 
                            c("CHR_ID", "CHR_POS", "CHR_POS", "SNPS")]

write.table(obesity_bed, 
            paste0(main_filepath, "GWASCatalog/gwascat_obesity_associations.bed"),
            sep = "\t", row.names = F, quote = F, col.names = F)
