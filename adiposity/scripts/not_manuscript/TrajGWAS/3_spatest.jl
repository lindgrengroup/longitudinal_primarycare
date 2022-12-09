# Author: Samvida S. Venkatesh
# Date: 01/06/2022

cd()
pwd()

using DataFrames, CSV
using Statistics
using TrajGWAS
using WiSER
using LinearAlgebra
using BGEN
using Serialization

# number of threads to use? try a couple of values
BLAS.set_num_threads(1)

RESDIR = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/results/"
BGENDIR = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/IMPUTATION/"
SAMPLEFILE = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample"

# Get arguments
STRATA = ARGS[1] 
CHR = ARGS[2] 
#CHUNKIDX = parse(Int, ARGS[3])
#NCHUNKS = parse(Int, ARGS[4])

fitted_null = string(RESDIR, STRATA, "/null_model.jls")
pvalfile = string(RESDIR, STRATA, "/trajgwas_sumstats_chr", CHR, ".txt")
#pvalfile = string(RESDIR, STRATA, "/trajgwas_sumstats_chr", CHR, "_chunk", CHUNKIDX, ".txt")

# Read data
nm = open(deserialize, fitted_null)
genetic_iids_subsample = nm.ids

bgenfilename = string(BGENDIR, "ukb_imp_chr", CHR, "_v3")
mfifilename = string(BGENDIR, "ukb_mfi_chr", CHR, "_v3.txt")

ukb_data = Bgen(bgenfilename * ".bgen"; sample_path = SAMPLEFILE)
genetic_iids = map(x -> parse(Int, split(x, " ")[1]), samples(ukb_data))

# Match up samples between null model and genetic data
order_dict = Dict{Int, Int}()
for (i, iid) in enumerate(genetic_iids)
    order_dict[iid] = i
end

sample_indicator = falses(length(genetic_iids))
for v in genetic_iids_subsample
    sample_indicator[order_dict[v]] = true # extract only the samples being used for the analysis
end

# GWAS for each chromosome

## pre-filtering SNPs not passing the criteria 
## (MAF > 0.01, info score > 0.8, HWE pval 1E-6)
min_maf = 0.01
min_info_score = 0.8
min_hwe_pval = 1e-6
mfi = CSV.read(mfifilename, DataFrame; header=false)
mfi.Column8 = map(x -> x == "NA" ? NaN : parse(Float64, x), mfi.Column8) # Column8: info score
snpmask = (mfi.Column6 .> min_maf) .& (mfi.Column8 .> min_info_score) # Column6: MAF

## compute range to run the analysis
#chunksize = n_variants(ukb_data) ÷ NCHUNKS + (n_variants(ukb_data) % NCHUNKS > 0 ? 1 : 0)
#startidx = chunksize * (CHUNKIDX - 1) + 1
#endidx = min(chunksize * CHUNKIDX, n_variants(ukb_data))
#snpmask = snpmask[startidx:endidx]

#println("running for variants $startidx to $endidx")

# rearrange data in nm so that it matches bgen data
nullinds = indexin(genetic_iids[sample_indicator], nm.ids)
nm.obswts .= isempty(nm.obswts) ? nm.obswts : nm.obswts[nullinds]
nm.ids .= nm.ids[nullinds]
nm.nis .= nm.nis[nullinds]
nm.data .= nm.data[nullinds]
@assert genetic_iids[sample_indicator] == nm.ids "there is some issue -- sampleids not matching"
    
trajgwas(nm, bgenfilename * ".bgen", count(sample_indicator);
    samplepath=SAMPLEFILE,
    pvalfile=pvalfile,
    snpinds=snpmask,
    min_hwe_pval = min_hwe_pval,
    bgenrowinds = sample_indicator,
    #startidx = startidx,
    #endidx = endidx,
    usespa=true)
