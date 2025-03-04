## load libraries
library(tidyverse)
library(Rcpp)
library(parallel)

## set seed 
set.seed(4578)

## create output directory
dir.create("data")

sourceCpp("TruncSkellams.cpp")

sourceCpp("tnorm.cpp")

## source function to run PF and return log-likelihood
source('BPF.R')

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave1/disease.rds")) %>%
  select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## select the 10th of the input parameters for simulation
inputs = readRDS("wave1/inputs.rds")
#You can select custom input parameters for the simulation by modifying 
#the parameter index value (idx_par) defined below
idx_par = 10 
saveRDS(inputs[idx_par,], "data/trueParams.rds")

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
N <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- N
u[2, ] <- I0

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- PF(pars[idx_par,], C = contact, u = u, ndays = 90, npart = 50, MD = TRUE, a_dis = 0.05, b_dis = 0.01, PF = FALSE)

## set stage names
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)

## collapse to data frame
disSims <- map(disSims[[1]], function(x, nms) {
        map(1:length(x), function(i, x, nms) {
            t(x[[i]]) %>%
            as.vector() %>%
            c(i)
        }, x = x, nms = nms) %>%
        {do.call("rbind", .)}
    })
disSims <- map(1:length(disSims), function(i, x) {
        cbind(x[[i]], i)
    }, x = disSims) %>%
    {do.call("rbind", .)}
colnames(disSims) <- c(stageNms, "rep", "t")
disSims <- as_tibble(disSims)

## extract simulation closest to median
medRep <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = (median - n)^2) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(1) %>%
    inner_join(disSims, by = "rep") %>%
    dplyr::select(!c(diff, rep))

## Skellam noise model for observations
a1 <- 0.01
a2 <- 0.2
b <- 0.1

skelNoise <- function(count, a1, a2, b1, b2) {
    ## create incidence over time
    inc <- diff(c(0, count))
    ## add observation noise
    for(i in 1:length(inc)) {
        inc[i] <- inc[i] +  rtskellam_cpp( a1 + b1 * inc[i], a2 + b2 * inc[i], -inc[i], 100000000)
    }
    ## return cumulative counts
    cumsum(inc)
}

skelNoise_H <- function(count, a1, a2, b1, b2) {
  ## add observation noise
  for(i in 1:length(count)) {
    count[i] <- count[i] + rtskellam_cpp(a1 + b1 * count[i], a2 + b2 * count[i], -count[i], 100000000)
  }
  ## return  counts
  count

}

medRep <- mutate(medRep, across(starts_with("H"), skelNoise_H, a1 = a1, a2 = a2, b1 = b, b2 = b, .names = "{.col}obs")) 

medRep <- mutate(medRep, across(starts_with("DI"), skelNoise, a1 = a1, a2 = a2, b1 = b, b2 = b, .names = "{.col}obs")) %>%
    mutate(across(starts_with("DH"), skelNoise, a1 = a1, a2 = a2, b1 = b, b2 = b, .names = "{.col}obs"))

## plot replicates
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(dplyr::select(medRep, !ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
            col = "red", linetype = "dashed"
        ) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(dplyr::select(medRep, t, ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("obs", "", var)) %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
            col = "blue", linetype = "dashed"
        ) +
        facet_grid(var ~ age, scales = "free") +
        xlab("Days") + 
        ylab("Counts")
ggsave("data/sims.pdf", p, width = 10, height = 10)

saveRDS(medRep, "data/disSims.rds")
saveRDS(pars, "data/pars.rds")

