###############################################################
## Plotting NROY trajectories code
###############################################################

### Setup ###

library(tidyverse)
library(patchwork)

library(tictoc)
source("inputs/dataTools.R")
sourceCpp("discreteStochModel.cpp")
sourceCpp("TruncSkellams.cpp")
sourceCpp("tnorm.cpp")
source("BPF_foreCasts.R")
library(parallel)


map <- purrr::map
### load contact matrix ###
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))
# smart_round function
smart_round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

###### Load data used to generate NROY spaces
## Note that if there is no data/disSims.rds, you will need to generate it by running source('stochModel.R') but this is not a good idea given that the data generated will not be the same as that used for NROY space
data <- readRDS("data/disSims.rds")
wave_idx = 10 # which wave (ensemble) you want to plot the trajectories
pf_days <- 48 #Ensure this is how many days of data were assimilated in the calibration
predict_days <- 30 #How far forward to go. I figured predicting 1 month ahead was interesting, but can be changed.
npoints <- 250
storage_directory= paste0("trajectoriesWave",wave_idx)
a_dis=0.05
b_dis=0.01 
num_particles=500

##### Preliminaries for the Particle Filter
### Code for inputs-disease files ###
# read in the disease file located in wavex folder
pars <- readRDS(paste0("wave", wave_idx, "/disease.rds")) %>%
  select(nu, nuA, !output)

## set up number of initial individuals in each age-class
N <- 10000
N <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 0)
S0 <- N - I0
## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- S0
u[2, ] <- I0

## Create a directory for particle data files to be stored, suggest using the name of the design file or similar
dir.create(storage_directory)
Nruns <- dim(pars)[1]
n_inbatch <- Nruns/10 #assume most designs are multiples of ten...

#Run the particle filter in batches to save memory.
batch <- 1
while(n_inbatch*batch <= Nruns){
  print(batch)
  tic()
  runs_md <- PF_sim(pars[(1+(batch-1)*n_inbatch):(batch*n_inbatch), ], C = C, data = data, u = u, ndays = pf_days, simdays = predict_days, npart = num_particles, MD = TRUE, a_dis = a_dis, b_dis = b_dis, saveAll = TRUE)

  # save counts data
  sims_md <- map_dfr(runs_md$particles, function(e) {map(e,  ~map(., ~as.vector(t(.)))) %>%
      map(~do.call("rbind", .)) %>%
      map(as_tibble) %>%
      bind_rows(.id = "t")}, .id= "run") %>%
    mutate(run= (batch-1)*n_inbatch + as.numeric(run))
  sims_md <- sims_md %>% mutate(t=as.numeric(t))
  stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
  colnames(sims_md) <- c("run", "t", stageNms)
  saveRDS(sims_md, file=paste(storage_directory, "/sims_batch_",batch,"counts.RDS", sep=""))
  rm(runs_md)
  rm(sims_md)
  toc()
  batch <- batch+1
}

######################################## 
## Plotting
######################################## 

#We can collate all runs for plotting first via
ensembleDesign <-  read_csv(paste0("wave", wave_idx,'/inputsWave', wave_idx, '.csv'))
ensembleDesign  = mutate(ensembleDesign , Plausible = rep(TRUE, dim(ensembleDesign)[1]))

sims_md <- list.files(path=storage_directory, pattern="counts.RDS",full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  mutate(NROY = ensembleDesign$Plausible[run]) %>%
  filter(NROY) 
#The following code plots the trajectories. The notPlot vector tells R which outputs to leave out. Sometimes if having memory issues, you have to remove a few of these.
notPlot <- c("RA", "RI", "RH", "P" ,"S")

age_label <- as.character(1:8)
age_label <- c("<5", "5-17", 
               "18-29", "30-39", "40-49", "50-59", 
               "60-69", "70+")
names(age_label) <- as.character(1:8)
predict_days <- 12 #How far forward to go. I figured predicting 1 month ahead was interesting, but can be changed.
p1 <-  dplyr::select(sims_md, -starts_with(notPlot)) %>%
  filter(t <= pf_days+predict_days)%>%
  pivot_longer(!c(run,t, NROY), names_to = "var", values_to = "n") %>%
  group_by(t, var, NROY) %>%
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
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = "50% Credible interval"), alpha = 0.5) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = "95% Credible interval"), alpha = 0.5) +
  geom_line(aes(y = median, colour = "Median")) +
  geom_line(
    aes(y = n, colour="Truth"),
    data = pivot_longer(select(data, !ends_with("obs")) %>% select(-starts_with(notPlot))%>% select(!ends_with("cum")) %>% filter(t <= pf_days+predict_days), !t, names_to = "var", values_to = "n") %>%
      mutate(age = gsub("[^0-9]", "", var)) %>%
      mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
      mutate(var = gsub("one", "1", var)) %>%
      mutate(var = gsub("two", "2", var)), linetype = "dashed") +
  geom_line(
    aes(y = n, colour = "Observations"),
    data = pivot_longer(select(data, t, starts_with(c("DH","DI")) &ends_with("obs")) %>% select(-starts_with(notPlot))%>% select(!ends_with("cumobs")) %>% filter(t <= pf_days+predict_days), !t, names_to = "var", values_to = "n") %>%
      mutate(var = gsub("obs", "", var)) %>%
      mutate(age = gsub("[^0-9]", "", var)) %>%
      mutate(var = gsub("[^a-zA-Z]", "", var)), linetype = "dashed") +
  geom_vline(aes(xintercept=pf_days, colour = "End of DA"), linetype="dotted") +
  theme(legend.position = "right") + scale_fill_manual("",values = c("#009ac4", "#00BFC4"))  + scale_colour_manual("", values = c("Median" = "#009ac4","Truth" = "red","Observations" ="blue", "End of DA" ="black" )) +
  facet_grid(var ~ age , scales = "free_y", labeller = labeller(age = age_label)) +
  xlab("Days") +
  ylab("Counts")+
  ggtitle(paste0("Wave", wave_idx, "trajectories"))
# p
ggsave('perfectModel_Ndays48_D_95CI.pdf', p1, width = 10)
