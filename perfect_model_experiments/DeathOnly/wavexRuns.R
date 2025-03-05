## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(tictoc)

## can run on the command line e.g.
## R CMD BATCH --no-save --no-restore --slave '--args 1' wavexRuns.R

## the first argument is the wave idx number

# check if being run in batch mode
args <- commandArgs(TRUE)
if(length(args) != 0) {
    ## extract command line arguments
    args <- commandArgs(TRUE)
    if(length(args) > 0) {
        stopifnot(length(args) == 1)
        wave_idx <- as.numeric(args[1])
    } else {
        stop("No arguments")
    }
} else {
    wave_idx <- 1
}
print(wave_idx)

## source simulation model
sourceCpp("discreteStochModel.cpp")

## read in truncated Skellam sampler
# source("trSkellam.R")
sourceCpp("TruncSkellams.cpp")

# source("RversionTNorms.R")
sourceCpp("tnorm.cpp")

## source function to run PF and return log-likelihood
source('BPF.R')

## set wave
wave <- wave_idx

## read in simulated data and generate incidence curves
data <- readRDS("data/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## set up number of initial individuals in each age-class
N <- 10000
N <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 0)
S0 <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- S0
u[2, ] <- I0

## run PF with some model discrepancy
tic()
runs_md <- PF(pars, C = contact, data = data, u = u, ndays = 48, npart = 500, MD = TRUE, 
    a_dis = 0.05, b_dis = 0.01, a1 = 0.01, a2 = 0.2, b = 0.1, saveAll = NA)
toc()

## save outputs
saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
