# Author: Samvida S. Venkatesh
# Date: 05/12/2022
# ADAPTED FROM: Frederik Heymann Lassen

library(argparse)
library(data.table)

main <- function(args){
  
  stopifnot(file.exists(args$in_file)) 
  d <- fread(args$in_file)
  genes <- fread('/well/lindgren/flassen/software/magma/auxiliary_files/genes/GRCh37/NCBI37.gene.loc')
  gene_mapping <- data.table(GENE = genes$V1, HGNC_SYMBOL=genes$V6)
  mrg <- merge(d, gene_mapping, all.x = TRUE)
  mrg$FDR <- stats::p.adjust(mrg$P, method = 'fdr')
  mrg <- mrg[order(mrg$P),]
  fwrite(mrg, file = args$out_file, sep = args$out_sep)
  
}


# add arguments
parser <- ArgumentParser()
parser$add_argument("--in_file", default=NULL, required = TRUE, help = "Input")
parser$add_argument("--out_file", default=NULL, required = TRUE, help = "Output")
parser$add_argument("--out_sep", default="\t", help = "Output seperator")
args <- parser$parse_args()

main(args)

