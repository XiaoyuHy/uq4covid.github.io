## load libraries
# X11()
# dev.off()
library(MASS)
library(tidyverse)
library(hmer)  ## experimental branch needed
library(mclust)
library(fields)
library(dgpsi) ## needs development version for set_seed()
init_py()

## can run on the command line e.g.
## R CMD BATCH --no-save --no-restore --slave '--args 2' emulateDGP_wavex.R
## the argument 2 is the new wave number

## check if being run in batch mode
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
}

## read in emulator objects
source("emulator_fns.R")

#########################################
###   NEED TO UPDATE VECTOR OF PATHS  ###
###   BELOW                           ###
#########################################
#wave_idx = 9
## set vector of paths to previous waves (including current wave)
prev_waves <- paste0("wave", 1:wave_idx)

#current wave folder
current_wave_folder <- paste0("wave", wave_idx)

## read in seeds
seeds <- readRDS(paste0(prev_waves[1], "/seeds.rds"))

## build custom emulator for previous waves
emulator <- list()
for(i in 1:(length(prev_waves) - 1)) {

    ## set seeds for reproducibility
    m_gp <- read(paste0(prev_waves[i], "/m_gp.pkl"))
    dgpsi::set_seed(seeds[i])
    m_gp <- set_imp(m_gp, B = 30)
    
    ## set up emulator
    emulator[[i]] <- list(ll = Proto_emulator$new(
            readRDS(paste0(prev_waves[1], "/ranges.rds")),
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
            range_list = readRDS(paste0(prev_waves[1], "/ranges.rds")), 
            epsilon = 0.00001,
            actives = readRDS(paste0(prev_waves[i], "/actives.rds")),
            target = readRDS(paste0(prev_waves[i], "/targets.rds")),
            pcutoff = readRDS(paste0(prev_waves[i], "/cutoff.rds")),
            voronoi = readRDS(paste0(prev_waves[i], "/voronoi.rds"))
        )
    )
}

## read in original training data across waves
inputs <- list()
ll <- list()
for(i in 1:(length(prev_waves) - 1)) {
    inputs[[i]] <- readRDS(paste0(prev_waves[i], "/inputs.rds"))
    ll[[i]] <- readRDS(paste0(prev_waves[i], "/runs_md.rds"))
}
inputs <- reduce(inputs, rbind)
ll <- reduce(ll, c)
inputs <- select(inputs, !output)

## check against implausibility at previous waves
## and extract non-implausible points across
## all waves as training points
new_set <- nth_implausible(emulator, inputs, list(ll = list(val = 1, sd = 1)), n = 1)
prev_points <- inputs[new_set <= 1, ]
prev_ll <- ll[new_set <= 1]

## check best point is in NROY
temp <- map_dbl(prev_waves[1:(length(prev_waves) - 1)], ~{
    readRDS(paste0(., "/targets.rds"))$ll$val
})
if(length(temp) > 1) {
    stopifnot(all(diff(temp) >= 0))
    temp <- max(temp)
}   
stopifnot(max(-exp(log(-prev_ll))) == temp)

## read in current runs
inputs <- readRDS(paste0(current_wave_folder, "/inputs.rds"))
ll <- readRDS(paste0(current_wave_folder, "/runs_md.rds"))

## extract training and test data
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
    
## bind to previous runs
if(length(prev_ll) > 0) {
    training <- rbind(training, as.matrix(cbind(prev_points, ll = log(-prev_ll))))
}
print(dim(training)[1])
## separate inputs and outputs
training_ll <- training[, ncol(training)]
training <- training[, -ncol(training)]
test_ll <- test[, ncol(test)]
test <- test[, -ncol(test)]

## scale inputs to (0, 1)
ranges <- readRDS(paste0(prev_waves[1], "/ranges.rds"))
for(i in 1:ncol(training)) {
    temp_ranges <- ranges[[colnames(training)[i]]]
    training[, i] <- (training[, i] - temp_ranges[1]) / diff(temp_ranges)
    test[, i] <- (test[, i] - temp_ranges[1]) / diff(temp_ranges)
}
stopifnot(all(apply(training, 2, function(x) all(x > 0 & x < 1))))
stopifnot(all(apply(test, 2, function(x) all(x > 0 & x < 1))))

## set target
targets <- list(ll = list(val = max(-exp(c(training_ll, test_ll))), sigma = 1))

#########################################
###   THE PART BELOW SOMETIMES NEEDS  ###
###   MANUAL INTERVENTION IF SAY      ###
###   SETTING SOME VARIABLES AS       ###
###   INACTIVE CAUSE OTHERS TO BECOME ###
###   INACTIVE CAUSE                  ###
#########################################

# combine scaled training and test data to perform sensivitivy analyis
tmp_ll = training_ll
training_all_temp = cbind(training, tmp_ll)
tmp_ll = test_ll
test_all_temp = cbind(test, tmp_ll)
allData = rbind(training_all_temp, test_all_temp)

#perform sensivitivy analyis
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
par(mar=c(7,4,4,2),mfrow=c(1,1))
figFileName = paste0(current_wave_folder, '/sensitivityWav', wave_idx, '.png', sep='')
png(file=figFileName,width=600, height=350)
barplot(mainEffects,names.arg = names(Sensitivity)[c(1:13)],las=2,main= paste('Wave ', wave_idx, sep=''))
dev.off()

## find the first 10 active variables 
names <- names(Sensitivity)[c(1:13)][order(mainEffects,decreasing = TRUE)]
main_effects = names[1:10]
colnames(training)
idx = match(main_effects, colnames(training))
actives = rep(FALSE,ncol(training))
actives[idx] = TRUE
## save active variables
saveRDS(actives,  paste0(current_wave_folder, "/actives.rds"))

## AT THIS POINT YOU SHOULD HAVE ACTIVE VARIABLES
colnames(training)[actives]

## save emulator
saveRDS(targets, paste0(current_wave_folder, "/targets.rds"))

## fit two-layer deep GP
m_gp <- dgp(
    training[, actives],
    training_ll, 
    name = "matern2.5",
    depth = 2,
    likelihood = "Hetero",
    N = 3000,
    B = 30
)
# X11()
# trace_plot(m_gp, layer = 1, node = 1)
# X11()
# trace_plot(m_gp, layer = 1, node = 2)
pdf(paste0(current_wave_folder, "/traces.pdf"))
trace_plot(m_gp, layer = 1, node = 1)
trace_plot(m_gp, layer = 1, node = 2)
dev.off()

## save object
write(m_gp, paste0(current_wave_folder, "/m_gp.pkl"))

#########################################
###   DGP USES LAST HALF OF CHAIN     ###
###   FOR INFERENCE, SO MUST HAVE     ###
###   CONVERGED BEFORE HALFWAY        ###
###   ELSE NEEDS RUNNING FOR LONGER   ###
###   (SEE BELOW)                     ###
#########################################

### if you need to run for longer, then 
### continue training
# m_gp <- continue(m_gp, N = 1500)
# X11()
# trace_plot(m_gp, layer = 1, node = 1)
# X11()
# trace_plot(m_gp, layer = 1, node = 2)
# pdf(paste0(current_wave_folder, "/traces.pdf"))
# trace_plot(m_gp, layer = 1, node = 1)
# trace_plot(m_gp, layer = 1, node = 2)
# dev.off()
# 
# ## save object
# write(m_gp, paste0(current_wave_folder, "/m_gp.pkl"))

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
pdf(paste0(current_wave_folder, "/validation.pdf"))
plot(m_gp)
plot(
    m_gp, 
    test[, actives],
    test_ll
)
dev.off()

## add test points into emulator object
m_gp <- update(
    m_gp, 
    X = rbind(training[, actives], test[, actives]), 
    Y = c(training_ll, test_ll), 
    B = 1
)

## save updated object
write(m_gp, paste0(current_wave_folder, "/m_gp.pkl"))

## set cutoff
saveRDS(0.95,  paste0(current_wave_folder, "/cutoff.rds"))

## CHECK FOR DOUBT POINTS AND SET VORONOI REGIONS IF REQUIRED

## set seeds for reproducibility
dgpsi::set_seed(seeds[length(prev_waves)])
m_gp <- set_imp(m_gp, B = 30)

## extract failed points
imp_mean <- predict(m_gp, m_gp$data$X)$results$mean
imp_var <- predict(m_gp, m_gp$data$X)$results$var
XD <- pnorm(log(-log(0.00001) - targets$ll$val), mean = imp_mean, sd = sqrt(imp_var))
XF <- which((1 - XD) > readRDS(paste0(current_wave_folder, "/cutoff.rds")))

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
saveRDS(voronoi, paste0(current_wave_folder, "/voronoi.rds"))

## build custom emulator for current wave
i <- length(prev_waves)

## set seeds for reproducibility
dgpsi::set_seed(seeds[i])
m_gp <- set_imp(m_gp, B = 30)

## set up emulator
emulator[[i]] <- list(ll = Proto_emulator$new(
        readRDS(paste0(prev_waves[1], "/ranges.rds")), 
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
        range_list = readRDS(paste0(prev_waves[1], "/ranges.rds")), 
        epsilon = 0.00001,
        actives = readRDS(paste0(prev_waves[i], "/actives.rds")),
        target = readRDS(paste0(prev_waves[i], "/targets.rds")),
        pcutoff = readRDS(paste0(prev_waves[i], "/cutoff.rds")),
        voronoi = readRDS(paste0(prev_waves[i], "/voronoi.rds"))
    )
)

## check best point is in NROY
temp <- list()
temp_ll <- list()
for(i in 1:length(prev_waves)) {
    temp[[i]] <- readRDS(paste0(prev_waves[i], "/inputs.rds"))
    temp_ll[[i]] <- readRDS(paste0(prev_waves[i], "/runs_md.rds"))
}
temp <- reduce(temp, rbind)
temp_ll <- reduce(temp_ll, c)
temp <- select(temp, !output)
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)
new_set <- nth_implausible(emulator, temp, list(ll = list(val = 1, sd = 1)), n = 1)
temp_ll <- temp_ll[new_set <= 1]
temp <- temp[new_set <= 1, ]
stopifnot(max(-exp(log(-temp_ll))) == targets$ll$val)

## read in baseline samples from previous non-implausible space
targets <- list(ll = list(val = 1, sd = 1))
inputs <- readRDS(paste0(prev_waves[length(prev_waves) - 1], "/baseline_samples.rds"))

## calculate space removed
sprem <- space_removal(emulator[[length(emulator)]], targets, inputs, cutoff = 1)
sprem
saveRDS(sprem, paste0(current_wave_folder, "/space_removed.rds"))

## extract previous design and update to extract which of these are still non-implausible
## with new emulator / targets
new_set <- nth_implausible(emulator, inputs, targets, n = 1)
inputs <- inputs[new_set <= 1, ]

## generate baseline set of points for the next design (and subsequent waves)
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = inputs)
saveRDS(baseline_points,  paste0(current_wave_folder, "/baseline_samples.rds"))

## generate new points using maximin
new_set <- maximin(baseline_points, 200, readRDS(paste0(prev_waves[1], "/ranges.rds")))
new_points <- baseline_points[new_set, ]
baseline_points <- baseline_points[-new_set, ]
new_set <- maximin(baseline_points, 50, readRDS(paste0(prev_waves[1], "/ranges.rds")))
new_points <- rbind(new_points, baseline_points[new_set, ])
rownames(new_points) <- 1:nrow(new_points)
baseline_points <- baseline_points[-new_set, ]

## write new runs
write_csv(new_points, paste0(current_wave_folder, "/inputsWave", wave_idx + 1, ".csv", sep=''))

## visualise new points
waves <- purrr::map(prev_waves, ~readRDS(paste0(., "/inputs.rds")))
waves[[length(waves) + 1]] <- data.frame(new_points)

ranges <- emulator[[1]]$ll$ranges
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

ggsave(paste0(current_wave_folder, '/new_design_wave', wave_idx + 1, '.pdf'), p, width = 18, height=10)


## generate baseline set of points for use in the next wave
baseline_points <- generate_new_design(emulator, 1000, targets, cutoff = 1, method = "slice", plausible_set = baseline_points)
saveRDS(baseline_points, paste0(current_wave_folder, "/baseline_samples.rds"))

