## load libraries
# X11()
# dev.off()
library(MASS)
library(tidyverse)
library(hmer)  
library(mclust)
library(fields)
library(dgpsi) 
init_py()

## read in emulator objects
source("emulator_fns.R")

## set parameter ranges
ranges <- list(
    nuA = c(0, 1),
    R0 = c(2, 4.5),
    TE = c(0.1, 2),
    TI1 = c(2.8, 4.5),
    TI2 = c(0.0001, 0.5),
    TP = c(1.2, 3),
    alphaTH = c(-1.4, 1.8),
    etaTH = c(0.001, 0.08),
    alphaEP = c(-5, 0),
    alphaI1D = c(-20, 0),
    alphaHD = c(-20, 0),
    alphaI1H = c(-5, 0),
    eta = c(0, 0.05)
)

## save ranges
saveRDS(ranges, "wave1/ranges.rds")

## test it out with log-likelihood estimation
inputs <- readRDS("wave1/inputs.rds")
ll <- readRDS("wave1/runs_md.rds")

## build BL emulators [NOTE: transformation is log(-log(likelihood))]
training <- mutate(inputs, ll = ll) %>%
    slice(1:200) %>%
    select(!output) %>%
    mutate(ll = log(-ll)) %>%
    as.matrix()

test <- mutate(inputs, ll = ll) %>%
    slice(-c(1:200)) %>%
    select(!output) %>%
    mutate(ll = log(-ll)) %>%
    as.matrix()

# separate inputs and outputs
training_ll <- training[, ncol(training)]
training <- training[, -ncol(training)]
test_ll <- test[, ncol(test)]
test <- test[, -ncol(test)]

## scale inputs to (0, 1)
for(i in 1:ncol(training)) {
    temp_ranges <- ranges[[colnames(training)[i]]]
    training[, i] <- (training[, i] - temp_ranges[1]) / diff(temp_ranges)
    test[, i] <- (test[, i] - temp_ranges[1]) / diff(temp_ranges)
}
stopifnot(all(apply(training, 2, function(x) all(x > 0 & x < 1))))
stopifnot(all(apply(test, 2, function(x) all(x > 0 & x < 1))))

## set target (NOTE: back transformation to log-likelihood)
targets <- list(ll = list(val = max(-exp(c(training_ll, test_ll))), sigma = 1))

#########################################
###   THE PART BELOW SOMETIMES NEEDS  ###
###   MANUAL INTERVENTION IF SAY      ###
###   SETTING SOME VARIABLES AS       ###
###   INACTIVE CAUSE OTHERS TO BECOME ###
###   INACTIVE CAUSE                  ###
#########################################

## combine scaled training and test data to perform sensivitivy analyis
tmp_ll = training_ll
training_all_temp = cbind(training, tmp_ll)
tmp_ll = test_ll
test_all_temp = cbind(test, tmp_ll)
allData = rbind(training_all_temp, test_all_temp)

##perform sensivitivy analyis
Sensitivity <- allData
Sensitivity %>% as_tibble %>% gather('parameter', 'value', -"tmp_ll") %>%
 ggplot(aes(x=value, y=`tmp_ll`)) +
 geom_point(alpha = 0.5) +
 facet_wrap(~parameter) +
 labs(y="llik", x='input') +
 geom_smooth(method = mgcv::gam,
             formula = y ~ s(x, bs = "tp"),
             fill = "red",
             method.args = list(method="GCV.Cp"))
Sensitivity =  as.data.frame(Sensitivity)
mainEffects <- c()
for(i in c(1:13)){
 gam1 <- mgcv::gam(Sensitivity$tmp_ll~s(as.matrix(Sensitivity)[,i]))
 mainEffects <- c(mainEffects,var(gam1$fitted)/var(Sensitivity$tmp_ll))
}
# plot the main effects and save the plot
wave_idx = 1
par(mar=c(7,4,4,2),mfrow=c(1,1))
figFileName = paste0('wave1/sensitivityWav', wave_idx, '.png', sep='')
png(file=figFileName,width=600, height=350)
barplot(mainEffects,names.arg = names(Sensitivity)[c(1:13)],las=2,main= paste('Wave ', wave_idx, sep=''))
dev.off()

## find the first 10 active variables 
names <- names(Sensitivity)[c(1:13)][order(mainEffects,decreasing = TRUE)]
main_10_effects = names[1:10]
idx = match(main_10_effects, colnames(training))
actives = rep(FALSE,ncol(training))
actives[idx] = TRUE
saveRDS(actives, "wave1/actives.rds")

## AT THIS POINT YOU SHOULD HAVE ACTIVE VARIABLES
colnames(training)[actives]

## fit two-layer deep GP with heteroscedastic likelihood
m_gp <- dgp(
    training[, actives],
    training_ll, 
    name = "matern2.5",
    depth = 2,
    likelihood = "Hetero",
    N = 5000,
    B = 30
)
## check trace plots
# X11()
# trace_plot(m_gp, layer = 1, node = 1)
# X11()
# trace_plot(m_gp, layer = 1, node = 2)
pdf("wave1/traces.pdf")
trace_plot(m_gp, layer = 1, node = 1)
trace_plot(m_gp, layer = 1, node = 2)
dev.off()

## save object
write(m_gp, "wave1/m_gp.pkl")

#########################################
###   DGP USES LAST HALF OF CHAIN     ###
###   FOR INFERENCE, SO MUST HAVE     ###
###   CONVERGED BEFORE HALFWAY        ###
###   ELSE NEEDS RUNNING FOR LONGER   ###
###   (SEE BELOW)                     ###
#########################################

### if you need to run for longer, then 
## continue training
#m_gp <- continue(m_gp, N = 1500)
#X11()
#trace_plot(m_gp, layer = 1, node = 1)
#X11()
#trace_plot(m_gp, layer = 1, node = 2)
#pdf("wave1/traces.pdf")
#trace_plot(m_gp, layer = 1, node = 1)
#trace_plot(m_gp, layer = 1, node = 2)
#dev.off()

## save object
#write(m_gp, "wave1/m_gp.pkl")

## check validation plots
summary(m_gp)
m_gp <- set_imp(m_gp, B = 30)
m_gp <- validate(m_gp)
m_gp <- validate(
    m_gp, 
    test[, actives],
    test_ll
)
# X11()
# plot(m_gp)
# X11()
# plot(
#     m_gp,
#     test[, actives],
#     test_ll
# )
pdf("wave1/validation.pdf")
plot(m_gp)
plot(
    m_gp,
    test[, actives],
    test_ll
)
dev.off()

## add test points into emulator
m_gp <- update(
    m_gp, 
    X = rbind(training[, actives], test[, actives]), 
    Y = c(training_ll, test_ll), 
    B = 1
)

## save updated object
write(m_gp, "wave1/m_gp.pkl")

## generate random seeds for reproducibility of
## imputations in future waves
seeds <- round(runif(30, 0, 50000000))
saveRDS(seeds, "wave1/seeds.rds")

## save variables for emulator
saveRDS(targets, "wave1/targets.rds")
saveRDS(0.95, "wave1/cutoff.rds")
# saveRDS(actives, "wave1/actives.rds")

## CHECK FOR DOUBT POINTS AND SET VORONOI REGIONS IF REQUIRED

## set seeds for reproducibility
dgpsi::set_seed(seeds[1])
m_gp <- set_imp(m_gp, B = 30)

## extract failed points
imp_mean <- predict(m_gp, m_gp$data$X)$results$mean
imp_var <- predict(m_gp, m_gp$data$X)$results$var
XD <- pnorm(log(-log(0.00001) - targets$ll$val), mean = imp_mean, sd = sqrt(imp_var))
XF <- which((1 - XD) > readRDS("wave1/cutoff.rds"))

## calculate doubt points
voronoi <- NULL
if(length(XF) > 0) {
    XD <- XF[-exp(m_gp$data$Y[XF, 1]) - targets$ll$val > log(0.00001)]

    if(length(XD) > 0) {
        ## extract lengthscale
        lscale <- m_gp$specs$layer1$node1$lengthscales

        ## augment space of doubt points if necessary
        temp <- matern(
            m_gp$data$X[-XD, , drop = FALSE],
            m_gp$data$X,
            2.5,
            lscale
        )
        temp <- which(apply(temp, 1, function(x, XD) {
            inds <- which(x == max(x, na.rm = TRUE))
            any(inds %in% XD)
        }, XD = XD))
        XD <- sort(unique(c(XD, temp)))
        if(length(XD) > 0) {
            voronoi <- list(kernel = matern, XD = XD)
        }
    }
}   

## save Voronoi region information
saveRDS(voronoi, "wave1/voronoi.rds")

## build custom emulator
emulator <- list(ll = Proto_emulator$new(
        readRDS("wave1/ranges.rds"),
        "ll",
        pred_func,
        pred_var_func,
        implausibility,
        print_func = function(em) {
            print(summary(em))
        },
        em = m_gp,
        hospStays = hospStays,
        hospThresh = hospThresh,
        pathwaysMod = pathwaysMod,
        pathwaysThresh = pathwaysThresh,
        ages = ages,
        pathwaysLimitFn = pathwaysLimitFn,
        NGM = NGM, 
        N = N, 
        S0 = S0, 
        C = C,
        range_list = readRDS("wave1/ranges.rds"),
        epsilon = 0.00001,
        actives = readRDS("wave1/actives.rds"),
        target = readRDS("wave1/targets.rds"),
        pcutoff = readRDS("wave1/cutoff.rds"),
        voronoi = readRDS("wave1/voronoi.rds")
    )
)

## check best point is in NROY
temp <- readRDS("wave1/inputs.rds")
temp_ll <- readRDS("wave1/runs_md.rds")
temp <- select(temp, !output)
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)
new_set <- nth_implausible(emulator, temp, list(ll = list(val = 1, sd = 1)), n = 1)
temp_ll <- temp_ll[new_set <= 1]
temp <- temp[new_set <= 1, ]
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)

## generate large sample from input space

## source dataTools
source("inputs/dataTools.R")
library(lhs)

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", "nuA"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1),
    stringsAsFactors = FALSE
) 

## generate LHS design
ndesign <- 1000
design <- randomLHS(ndesign, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
inputs <- convertDesignToInput(design, parRanges, "zero_one")

## generate space-filling design for other parameters

## generate design points for hospital stay lengths
hospStaysInput <- FMMmaximin(
    hospStays, 
    ndesign + 20, 
    matrix(c(-Inf, Inf, 0, Inf), ncol = 2, byrow = TRUE)
) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)
    
## check against prior density restrictions
hospStaysInput <- hospStaysInput[dens(as.matrix(hospStaysInput), hospStays$modelName, hospStays$parameters, logarithm = TRUE) > hospThresh, ]
if(nrow(hospStaysInput) < ndesign) stop("Can't generate enough valid hospital points")
hospStaysInput <- hospStaysInput[1:ndesign, ]

## generate design points for other transition probabilities

## produces design points subject to constraints
pathwaysInput <- FMMmaximin(
    pathwaysMod, 
    ndesign + 20,
    matrix(c(rep(c(-20, 0), times = 4), 0, 1), ncol = 2, byrow = TRUE),
    pathwaysLimitFn,
    ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)
) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5)
    
## check against prior density restrictions
pathwaysInput <- pathwaysInput[dens(as.matrix(pathwaysInput), pathwaysMod$modelName, pathwaysMod$parameters, logarithm = TRUE) > pathwaysThresh, ]
if(nrow(pathwaysInput) < ndesign) stop("Can't generate enough valid pathways points")
pathwaysInput <- pathwaysInput[1:ndesign, ]

## bind to design
inputs <- cbind(inputs, hospStaysInput, pathwaysInput)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = "Ens1", nrow(inputs))

## convert input to disease
disease <- convertInputToDisease(inputs, C, N, S0, ages)
inputs <- semi_join(inputs, disease, by = "output")
    
## bind to design
inputs <- select(inputs, !!colnames(training))
saveRDS(inputs, "wave1/initial_samples.rds")

## calculate space removed
targets <- list(ll = list(val = 1, sd = 1))
sprem <- space_removal(emulator, targets, inputs, cutoff = 1)
sprem
saveRDS(sprem, "wave1/space_removed.rds")

## extract previous design and update to extract which of these are still non-implausible
## with new emulator / targets
new_set <- nth_implausible(emulator, inputs, targets, n = 1)
inputs <- inputs[new_set <= 1, ]

## generate baseline set of points for the next design (and subsequent waves)
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = inputs)
saveRDS(baseline_points, "wave1/baseline_samples.rds")

## generate new points using maximin
new_set <- maximin(baseline_points, 200, readRDS("wave1/ranges.rds"))
new_points <- baseline_points[new_set, ]
baseline_points <- baseline_points[-new_set, ]
new_set <- maximin(baseline_points, 50, readRDS("wave1/ranges.rds"))
new_points <- rbind(new_points, baseline_points[new_set, ])
rownames(new_points) <- 1:nrow(new_points)
baseline_points <- baseline_points[-new_set, ]

## write new runs
write_csv(new_points, "wave1/inputsWave2.csv")

## visualise new points
waves <- list(readRDS("wave1/inputs.rds"))
waves[[length(waves) + 1]] <- new_points

ranges <- emulator$ll$ranges
ranges <- ranges[names(ranges) %in% colnames(waves[[1]])]

p <- wave_points(waves, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
                 wave_numbers = 1:length(waves)) +
  labs(fill = "Wave", title = NULL)
## set path to outputs
pathtofolder <- "./"
#add true parameters to the plot
pars1 <- readRDS(paste0(pathtofolder, "outputs/trueParams.rds")) 
## extract scales
g <- ggplot_build(p[1, 1])
fill_cols <- c(unique(g$data[[1]]["fill"])$fill, "red")
pt_cols <- c(unique(g$data[[1]]["fill"])$fill, "red")
waves[[length(waves) + 1]] <- pars1
p <- wave_points(waves, input_names = names(ranges), p_size = 1, zero_in = FALSE, 
                 wave_numbers = 1:length(waves)) +
  labs(fill = "Wave", title = NULL) +
  scale_fill_manual(values = fill_cols, labels = c(1:(length(waves) - 1), "Truth")) +
  scale_colour_manual(values = pt_cols)            

ggsave(paste0('wave1/new_design_wave', wave_idx+1, '.pdf'), p, width = 18, height = 10)

## generate baseline set of points for use in the next wave
## THIS CAN TAKE A LONG TIME AT LATER WAVES BUT YOU CAN START
## THE NEW WAVE OFF TO RUN AND THEN RETURN TO THIS - THIS PROVIDES
## SAMPLES USED TO EVALUATE THE CUT-OUT SPACE AND AS A BASIS FOR
## THE SLICE SAMPLER IN LATER WAVES
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = baseline_points)
saveRDS(baseline_points, "wave1/baseline_samples.rds")

