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
#data <- readRDS("data/disSims.rds") # This line is commented out since the real data is loaded by lines 42-54.
wave_idx = 10 # which wave (ensemble) you want to plot the trajectories
pf_days <- 48 #Ensure this is how many days of data were assimilated in the calibration
predict_days <- 30 #How far forward to go. I figured predicting 1 month ahead was interesting, but can be changed.
npoints <- 250
storage_directory= paste0("trajectoriesWave",wave_idx)
a_dis=0.05
b_dis=0.01 
num_particles=500

## read in real data 
dataD <- readRDS("deathTotalByAge.rds")
# total hosptialisation cases
dataH <- readRDS("hospitalCasesTotal_Eng_forDemo.rds")
dataH <- select(dataH, -date)

#Cumulative deaths observations
cumDobs <- dataD 
for (i in (2:length(dataD$t))){
  cumDobs[i,1:8] <- colSums(dataD[1:i,1:8]) 
}
data <- list(dataD,dataH)
names(data) = c("dataD", "dataH")

##### Preliminaries for the Particle Filter
### Code for inputs-disease files ###
# read in the disease file located in wavex folder
pars <- readRDS(paste0("wave", wave_idx, "/disease.rds")) %>%
  select(nu, nuA, !output)

## set up number of initial individuals in each age-class
N <- 56550138
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
  tic()
  runs_md <- PF_sim(pars[(1+(batch-1)*n_inbatch):(batch*n_inbatch), ], C = C, data = data, u = u, ndays = pf_days, simdays = predict_days, npart = num_particles, MD = TRUE, a_dis = a_dis, b_dis = b_dis, saveAll = TRUE)
  # save counts data
  sims_md <- map_dfr(runs_md$particles, function(e) {map(map(e,1),  ~map(., ~as.vector(t(.)))) %>%
      map(~do.call("rbind", .)) %>%
      map(as_tibble) %>%
      bind_rows(.id = "t")}, .id= "run") %>%
    mutate(run= (batch-1)*n_inbatch + as.numeric(run))
  sims_md <- sims_md %>% mutate(t=as.numeric(t))
  # print(sims_md$t)
  stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
  colnames(sims_md) <- c("run", "t", stageNms)
  saveRDS(sims_md, file=paste(storage_directory, "/sims_batch_",batch,"counts.RDS", sep=""))
  
  sims_mdHpre <- map_dfr(runs_md$particles, function(e) {map(map(e,2),  ~map(., ~as.vector(t(.)))) %>%
      map(~do.call("rbind", .)) %>%
      map(as_tibble) %>%
      bind_rows(.id = "t")}, .id= "run") %>%
    mutate(run= (batch-1)*n_inbatch + as.numeric(run))
  sims_mdHpre <- sims_mdHpre %>% mutate(t=as.numeric(t))
  
  stageNms <- c("Hpre")
  colnames(sims_mdHpre) <- c("run", "t", stageNms)
  saveRDS(sims_mdHpre, file=paste(storage_directory, "/sims_batch_",batch,"Hpre.RDS", sep=""))
  
  sims_mdCumD <- map_dfr(runs_md$particles, function(e) {map(map(e,3),  ~map(., ~as.vector(t(.)))) %>%
      map(~do.call("rbind", .)) %>%
      map(as_tibble) %>%
      bind_rows(.id = "t")}, .id= "run") %>%
    mutate(run= (batch-1)*n_inbatch + as.numeric(run))
  sims_mdCumD <- sims_mdCumD %>% mutate(t=as.numeric(t))
  
  stageNms <- map(c("CumD"), ~paste0(., 1:8)) %>%
    reduce(c)
  colnames(sims_mdCumD) <- c("run", "t", stageNms)
  saveRDS(sims_mdCumD, file=paste(storage_directory, "/sims_batch_",batch,"cumD.RDS", sep=""))
  rm(runs_md)
  rm(sims_md)
  rm(sims_mdHpre)
  rm(sims_mdCumD)
  toc()
  batch <- batch+1
}

######################################## 
## Plotting
######################################## 

#We can collate all runs for plotting first via
dataD <- readRDS("deathTotalByAge.rds")
ensembleDesign <-  read_csv(paste0("wave", wave_idx,'/inputsWave', wave_idx, '.csv'))
ensembleDesign  = mutate(ensembleDesign , Plausible = rep(TRUE, dim(ensembleDesign)[1]))
sims_md <- list.files(path=storage_directory, pattern="counts.RDS",full.names = TRUE) %>%
  map_dfr(readRDS) %>%
  mutate(NROY = ensembleDesign$Plausible[run]) %>%
  filter(NROY) 
#The following code plots the trajectories. The notPlot vector tells R which outputs to leave out. Sometimes if having memory issues, you have to remove a few of these.
notPlot <- c("RA", "RI", "RH", "P" ,"S")
p55 <-  dplyr::select(sims_md, -starts_with(notPlot))
p66 = mutate(p55,"D1" = rowSums(select(p55,c("DH1", "DI1"))), "D2" = rowSums(select(p55,c("DH2", "DI2"))),
             "D3" = rowSums(select(p55,c("DH3", "DI3"))), "D4" = rowSums(select(p55,c("DH4", "DI4"))), "D5" = rowSums(select(p55,c("DH5", "DI5"))),
             "D6" = rowSums(select(p55,c("DH6", "DI6"))), "D7" = rowSums(select(p55,c("DH7", "DI7"))),
             "D8" = rowSums(select(p55,c("DH8", "DI8")))) %>%select(!starts_with(c("DH", "DI")))

age_label <- as.character(1:8)
age_label <- c("<5", "5-17", 
               "18-29", "30-39", "40-49", "50-59", 
               "60-69", "70+")
names(age_label) <- as.character(1:8)
predict_days <- 14 #How far forward to go. I figured predicting 1 month ahead was interesting, but can be changed.
p1 <-  p66 %>%
  filter(t <= pf_days+predict_days)%>%
  pivot_longer(!c(run,t, NROY), names_to = "var", values_to = "n") %>%
  group_by(t, var, NROY) %>%
  summarise(
    LCI = quantile(n, probs = 0.05),
    LQ = quantile(n, probs = 0.25),
    median = median(n),
    UQ = quantile(n, probs = 0.75),
    UCI = quantile(n, probs = 0.95),
    .groups = "drop"
  ) %>%
  mutate(age = gsub("[^0-9]", "", var)) %>%
  mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
  mutate(var = gsub("one", "1", var)) %>%
  mutate(var = gsub("two", "2", var)) %>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = "50% Credible interval"), alpha = 0.5) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = "90% Credible interval"), alpha = 0.5) +
  geom_line(aes(y = median, colour = "Median")) +
  geom_line(
    aes(y = n, colour = "Observations"),
    data = pivot_longer(select(dataD, t, starts_with(c("D")) &ends_with("obs")) %>% select(-starts_with(notPlot))%>% select(!ends_with("cumobs")) %>% filter(t <= pf_days+predict_days), !t, names_to = "var", values_to = "n") %>%
      mutate(var = gsub("obs", "", var)) %>%
      mutate(age = gsub("[^0-9]", "", var)) %>%
      mutate(var = gsub("[^a-zA-Z]", "", var)), linetype = "dashed") +
  geom_vline(aes(xintercept=pf_days, colour = "Lockdown"), linetype="dotted") +
  theme(legend.position = "right") + scale_fill_manual("",values = c("#009ac4", "#00BFC4"))  + scale_colour_manual("", values = c("Median" = "#009ac4","Observations" ="blue", "Lockdown" ="black" )) +
  facet_grid(var ~ age , scales = "free_y", labeller = labeller(age = age_label)) +
  xlab("Days") +
  ylab("Counts")+
  ggtitle(paste0("Wave", wave_idx, "trajectories"))
# p
ggsave('realData_D_H_Betas1.pdf', p1, width = 10)

## Below code is for plotting trajectories of cumulative deaths and hospitalisations

## zoom solution from https://stackoverflow.com/questions/63550588/ggplot2coord-cartesian-on-facets
UniquePanelCoords <- ggplot2::ggproto(
  "UniquePanelCoords", ggplot2::CoordCartesian,
  
  num_of_panels = 1,
  panel_counter = 1,
  layout = NULL,
  
  setup_layout = function(self, layout, params) {
    self$num_of_panels <- length(unique(layout$PANEL))
    self$panel_counter <- 1
    self$layout <- layout # store for later
    layout
  },
  
  setup_panel_params =  function(self, scale_x, scale_y, params = list()) {
    train_cartesian <- function(scale, limits, name, given_range = c(NA, NA)) {
      if (anyNA(given_range)) {
        expansion <- ggplot2:::default_expansion(scale, expand = self$expand)
        range <- ggplot2:::expand_limits_scale(scale, expansion, coord_limits = limits)
        isna <- is.na(given_range)
        given_range[isna] <- range[isna]
      }
      out <- list(
        ggplot2:::view_scale_primary(scale, limits, given_range),
        sec = ggplot2:::view_scale_secondary(scale, limits, given_range),
        arrange = scale$axis_order(),
        range = given_range
      )
      names(out) <- c(name, paste0(name, ".", names(out)[-1]))
      out
    }
    
    
    this_layout <- self$layout[ self$panel_counter,, drop = FALSE ]
    self$panel_counter <- 
      if (self$panel_counter < self$num_of_panels) {
        self$panel_counter + 1
      } else 1
    
    # determine merge column names by removing all "standard" names
    layout_names <- setdiff(names(this_layout),
                            c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y"))
    limits_names <- setdiff(names(self$panel_limits),
                            c("xmin", "xmax", "ymin", "ymax"))
    
    limit_extras <- setdiff(limits_names, layout_names)
    if (length(limit_extras) > 0) {
      stop("facet names in 'panel_limits' not found in 'layout': ",
           paste(sQuote(limit_extras), collapse = ","))
    } else if (length(limits_names) == 0 && NROW(self$panel_limits) == 1) {
      # no panels in 'panel_limits'
      this_panel_limits <- cbind(this_layout, self$panel_limits)
    } else {
      this_panel_limits <- merge(this_layout, self$panel_limits, all.x = TRUE, by = limits_names)
    }
    
    if (isTRUE(NROW(this_panel_limits) > 1)) {
      stop("multiple matches for current panel in 'panel_limits'")
    }
    
    # add missing min/max columns, default to "no override" (NA)
    this_panel_limits[, setdiff(c("xmin", "xmax", "ymin", "ymax"),
                                names(this_panel_limits)) ] <- NA
    
    c(train_cartesian(scale_x, self$limits$x, "x",
                      unlist(this_panel_limits[, c("xmin", "xmax"), drop = TRUE])),
      train_cartesian(scale_y, self$limits$y, "y",
                      unlist(this_panel_limits[, c("ymin", "ymax"), drop = TRUE])))
  }
)

coord_cartesian_panels <- function(panel_limits, expand = TRUE, default = FALSE, clip = "on") {
  ggplot2::ggproto(NULL, UniquePanelCoords,
                   panel_limits = panel_limits,
                   expand = expand, default = default, clip = clip)
}


# sum the death and hospitalisations across the 8 ages
dataD[69:dim(dataD)[1],] <- NA
dataD <- mutate(dataD, 'Dobs' = rowSums(select(dataD, starts_with('D')))) %>%
  select(c(Dobs, t))
dataH[69:dim(dataH)[1],] <- NA
dataPlotCounts <-dataH[1:dim(dataD)[1],]
dataPlotInc <- dataD

cumDobs <- mutate(cumDobs, 'Dobs' = rowSums(select(cumDobs, starts_with('D')))) %>%
  select(c(Dobs, t))
cumDobs[69:dim(cumDobs)[1],] <- NA
ensembleDesign <-  read_csv(paste0("wave", wave_idx,'/inputsWave', wave_idx, '.csv'))
ensembleDesign  = mutate(ensembleDesign, Plausible = rep(TRUE, dim(ensembleDesign)[1]))
sims_mdHpre <- list.files(path="trajectories", pattern="Hpre.RDS",full.names = TRUE) %>%
  map_dfr(readRDS)
sims_md <- sims_mdHpre %>%
  mutate(NROY = ensembleDesign$Plausible[sims_mdHpre$run]) %>%
  filter(NROY) 

tmp2 <- sims_md %>%rename(H = "Hpre") %>% pivot_longer(!c(run,t, NROY), names_to = "var", values_to = "n") %>%
  mutate(var = gsub("[^a-zA-Z1]", "", var)) 
tmp3<- tmp2 %>% 
  group_by(t, var, NROY) %>%
  summarise(
    LCI = quantile(n, probs = 0.05),
    LQ = quantile(n, probs = 0.25),
    median = median(n),
    UQ = quantile(n, probs = 0.75),
    UCI = quantile(n, probs = 0.95),
    .groups = "drop"
  )
predict_days = 14
p1 <- tmp3 %>%
  filter(t <= pf_days+predict_days)%>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = "50% Credible interval"), alpha = 0.5) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = "90% Credible interval"), alpha = 0.5) +
  geom_line(aes(y = median, colour = "Median")) +
  geom_line(
    aes(y = n, colour = "Data"),
    data = pivot_longer(select(dataH, t, ends_with("obs"))%>% filter(t <= pf_days+predict_days), !t, names_to = "var", values_to = "n") %>%  mutate(var = gsub("obs", "", var)), linetype = "dashed"
  ) +
  geom_vline(aes(xintercept=pf_days, colour = "Lockdown"), linetype="dotted") +
  theme(legend.position = "right") + scale_fill_manual("",values = c("#009ac4", "#00BFC4"))  + scale_colour_manual("", values = c("Median" = "#009ac4","Data" ="blue", "Lockdown" ="red" )) +
  facet_grid(rows = vars(var) , scales = "free_y") +
  xlab("Days") +
  ylab("Hospitalisations")
# p1 

p2 <- p1 +
  coord_cartesian_panels(
    panel_limits = tibble::tribble(
      ~var,  ~ymin, ~ymax
      # , "D"           ,     0,     200
      , "H"           ,     0,     4000
      
    )
  )
# p2

p3 <- p1 + p2

sims_mdCumD <- list.files(path=storage_directory, pattern="cumD.RDS",full.names = TRUE) %>%
  map_dfr(readRDS)
sims_md <- sims_mdCumD %>%
  mutate(NROY = wave10Design$Plausible[sims_mdCumD$run]) %>%
  filter(NROY) 

tmp1 <- mutate(sims_md,'D' = rowSums(select(sims_md, starts_with('CumD')))) %>% select(c(run,t,NROY,D))
tmp2 <- tmp1 %>% pivot_longer(!c(run,t, NROY), names_to = "var", values_to = "n") %>%
  mutate(var = gsub("[^a-zA-Z1]", "", var)) 
tmp3<- tmp2 %>% 
  group_by(t, var, NROY) %>%
  summarise(
    LCI = quantile(n, probs = 0.05),
    LQ = quantile(n, probs = 0.25),
    median = median(n),
    UQ = quantile(n, probs = 0.75),
    UCI = quantile(n, probs = 0.95),
    .groups = "drop"
  )

p4 <- tmp3 %>%
  filter(t <= pf_days+predict_days)%>%
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = "50% Credible interval"), alpha = 0.5) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = "90% Credible interval"), alpha = 0.5) +
  geom_line(aes(y = median, colour = "Median")) +
  geom_line(
    aes(y = n, colour = "Data"),
    data = pivot_longer(select(cumDobs, t, ends_with("obs"))%>% filter(t <= pf_days+predict_days), !t, names_to = "var", values_to = "n") %>%  mutate(var = gsub("obs", "", var)), linetype = "dashed"
  ) +
  geom_vline(aes(xintercept=pf_days, colour = "Lockdown"), linetype="dotted") +
  theme(legend.position = "right") + scale_fill_manual("",values = c("#009ac4", "#00BFC4"))  + scale_colour_manual("", values = c("Median" = "#009ac4","Data" ="blue", "Lockdown" ="red" )) +
  facet_grid(rows = vars(var) , scales = "free_y") +
  xlab("Days") +
  ylab("Cumulative deaths")
# p1 

p5 <- p4 +
  coord_cartesian_panels(
    panel_limits = tibble::tribble(
      ~var,  ~ymin, ~ymax
      , "D"           ,     0,     1000
      # , "I1"           ,     0,     1000
      
    )
  )
# p2

p4 <- p4 +  ggtitle(paste0("Wave", wave_idx, "trajectories"))

p5 <- p5 + ggtitle(paste0("Zoomed-in view of wave", wave_idx, "trajectories"))

p6 <- p4 + p5

pp <- p6 /p3  + plot_layout(guides = 'collect')
ggsave('realData_D_H_cumulativeCounts_Betas1.pdf', pp, width = 12, height=6)


