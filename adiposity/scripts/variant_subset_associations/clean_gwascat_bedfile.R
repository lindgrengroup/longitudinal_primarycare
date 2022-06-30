# Author: Samvida S. Venkatesh
# Date: 29/06/22

library(tidyverse)
library(biomaRt)

main_filepath <- "/well/lindgren/samvida/Resources/GWASCatalog/"

bm <- useMart(biomart = "ENSEMBL_MART_SNP", 
              dataset = "hsapiens_snp")

# Read obesity associated SNPs  ----

obesity_snps <- read.table(paste0(main_filepath, "tmp_gwascat_obesity_associations.bed"),
                            sep = "\t", header = F, quote = "", stringsAsFactors = F)
colnames(obesity_snps) <- c("CHR", "POS0", "POS1", "SNP")

# Some SNP fields have multiple SNPs separated by "; " or " x "
all_snps <- unlist(strsplit(obesity_snps$SNP, "; "))
all_snps <- unlist(strsplit(all_snps, " x "))
obesity_snps <- obesity_snps %>% filter(!grepl(";|[ x ]", SNP))

not_in_dat <- all_snps[which(!all_snps %in% obesity_snps$SNP)]
add_in_dat <- data.frame(CHR = "", 
                         POS0 = as.numeric(""), 
                         POS1 = as.numeric(""), SNP = not_in_dat)
obesity_snps <- bind_rows(obesity_snps, add_in_dat)

# Get unique lines
obesity_snps <- obesity_snps %>% distinct()

# If missing chr/pos, add this in from SNP name
missing_dat <- obesity_snps %>% 
  filter(CHR == "" & grepl("^chr", SNP)) %>%
  mutate(tmp = gsub("chr", "", SNP),
         CHR = as.character(gsub(":.*", "", tmp)),
         POS0 = as.numeric(gsub(".*:", "", tmp)), 
         POS1 = POS0) 
missing_dat <- missing_dat[, c("CHR", "POS0", "POS1", "SNP")]

# If missing chr/pos but is an rsid, get this data from ENSEMBL
for_ensembl_check <- obesity_snps %>% 
  filter(CHR == "" & grepl("^rs", SNP))

ensembl_res <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
             filters = "snp_filter",
             values = for_ensembl_check$SNP,
             mart = bm)

for_ensembl_check <- for_ensembl_check %>%
  mutate(CHR = as.character(ensembl_res$chr_name[match(SNP, ensembl_res$refsnp_id)]),
         POS0 = as.numeric(ensembl_res$chrom_start[match(SNP, ensembl_res$refsnp_id)]),
         POS1 = POS0)

# Combine all the new mapped chr/pos files together
cleaned_bed <- obesity_snps %>%
  filter(CHR != "") %>%
  mutate(CHR = as.character(CHR),
         POS0 = as.numeric(POS0), POS1 = as.numeric(POS1))

cleaned_bed <- bind_rows(cleaned_bed, missing_dat)
cleaned_bed <- bind_rows(cleaned_bed, for_ensembl_check) 
cleaned_bed <- cleaned_bed[complete.cases(cleaned_bed), ]

# For bed formatting, change POS0 to POS1 - 1
cleaned_bed <- cleaned_bed %>%
  mutate(CHR = paste0("chr", CHR), 
         POS0 = POS1 - 1)

write.table(cleaned_bed, 
            paste0(main_filepath, "gwascat_obesity_associations_hg38.bed"),
            sep = "\t", row.names = F, quote = F, col.names = F)
