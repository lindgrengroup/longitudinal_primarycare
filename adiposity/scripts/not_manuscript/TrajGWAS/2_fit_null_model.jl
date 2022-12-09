# Author: Samvida S. Venkatesh
# Date: 31/05/2022

cd()
pwd()

# Load packages
using DataFrames, CSV
using Statistics
using TrajGWAS
using Ipopt, WiSER
using LinearAlgebra
using BGEN

# Get arguments
STRATA = ARGS[1]
INPUTDIR = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/hormone_ehr/TrajGWAS/data/"
OUTPUTDIR = "/well/lindgren-ukbb/projects/ukbb-11867/samvida/adiposity/TrajGWAS/results/"

# To fit the null model using a single thread
BLAS.set_num_threads(1)
solver = Ipopt.IpoptSolver(print_level=1, watchdog_shortened_iter_trigger=5, max_iter=120)

# Read data

dat_for_gwas = CSV.read(string(INPUTDIR, STRATA, "_pheno_covars.csv"), DataFrame)

path_nullresults = string(OUTPUTDIR, STRATA, "/null_results.txt")
path_nullmodel = string(OUTPUTDIR, STRATA, "/null_model.jls")

# Fit null model
# Fixed effects: intercept, age- and age-squared, 
# confounder adjusted - sex, genetic PCs
# Random effects: intercept and linear age term (slope)

if contains(STRATA, "sex_comb")
	@time nm = trajgwas(@formula(value ~ 1 + sex + age_event + age_sq + data_provider + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21),
    @formula(value ~ 1 + age_event),
    @formula(value ~ 1 + sex + age_event + age_sq + data_provider + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21),
    :eid, # subject ID
    dat_for_gwas,
    nothing;
    nullfile=path_nullresults,
    solver=solver,
    runs=10)
else
	@time nm = trajgwas(@formula(value ~ 1 + age_event + age_sq + data_provider + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21),
	@formula(value ~ 1),
    @formula(value ~ 1 + age_event + age_sq + data_provider + genotyping_array + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20 + PC21),
    :eid, # subject ID
    dat_for_gwas,
    nothing;
    nullfile=path_nullresults,
    solver=solver,
    runs=10)
end

println(nm)
using Serialization
open(path_nullmodel, "w") do io
    Serialization.serialize(io, nm)
end
