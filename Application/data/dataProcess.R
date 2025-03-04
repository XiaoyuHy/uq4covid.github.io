library(ggplot2)
library(tidyverse)
# load calibration data
# deaths by region and age band
# an array of dates x regions x PHE ages bands
# fl0 <- "https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv&release=2021-05-10"
dat0 <- read_csv('nation_2021-05-10.csv')

# ages we want (to avoid overlapping age bands)
ages <- c("00_04", "05_09", "10_14", "15_19", "20_24", "25_29", "30_34", 
          "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", 
          "70_74", "75_79", "80_84", "85_89", "90+")
dat1 <- subset(dat0, age %in% ages)
# dates
dates <- sort(unique(dat1$date))

# set as factors, to make life easy for putting into an array
dat1$age <- factor(dat1$age, levels = ages)
dat1$date <- as.factor(dat1$date)

# then put into array in form dates x regions x ages
# (maybe not the most logical order: use aperm to change)
dimnames <- list(date = dates,  age = ages)
arr <- array(0L, dim = sapply(dimnames, length), dimnames = dimnames)
ind <- cbind(as.integer(dat1$date), as.integer(dat1$age))
arr[ind] <- dat1$deaths
PHE_days <- dimnames(arr)[[1]]
n_PHE_days <- length(PHE_days)

# aggregate calibration data to MetaWards output age bands
n_ages <- 8
PHE_ages <- c("00_04", "05_09", "10_14", "15_19", "20_24", "25_29", "30_34", 
              "35_39", "40_44", "45_49", "50_54", "55_59", "60_64", "65_69", 
              "70_74", "75_79", "80_84", "85_89", "90+")
n_PHE_ages <- length(PHE_ages)
wts <- matrix(0, n_PHE_ages, n_ages)
wts[1, 1] <- 1 # 0-4 -> 0-4
wts[2:4, 2] <- c(1, 1, 3/5) # (5-9, 10-14, 15-19) -> 5-17
wts[4:6, 3] <- c(2/5, 1, 1) # (15-19, 20-24, 25-29) -> 18-29
wts[7:8, 4] <- c(1, 1) # (30-34, 35-39) -> 30-39
wts[9:10, 5] <- c(1, 1) # (40-44, 45-49) -> 40-49
wts[11:12, 6] <- c(1, 1) # (50-54, 55-59) -> 50-59
wts[13:14, 7] <- c(1, 1) # (60-64, 65-69) -> 60-69
wts[15:19, 8] <- c(1, 1, 1, 1, 1) # (70-74, 75-79, 80-84, 85-89, 90+) -> 70+
# then deaths
prep <- as.Date("2020-01-15") + 1:46
death <- array(NA, c(n_PHE_days, n_ages))
for (i in 1:n_PHE_days) {
  death[i, ] <- death0[i , ] %*% wts
}
death <- apply(death, 2, function(x) c(integer(length(prep)), x))

deathTotalByAge <- death
deathTotalByAge <- as.data.frame(deathTotalByAge) %>% mutate(t=c(2:477)) %>%
  rename(D1obs='V1', D2obs = 'V2', D3obs = 'V3', D4obs= 'V4', D5obs= 'V5', D6obs= 'V6', D7obs = 'V7', D8obs = 'V8')
saveRDS(deathTotalByAge, 'deathTotalByAge.rds')
