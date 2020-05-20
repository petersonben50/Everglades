#### code/geochemistry/geochem_gradients_2018.R ####
# Benjamin D. Peterson

# This script will generate the graphs I 
# want for highlighting the sulfur gradients
# along which we sampled in 2018.



#### Always start with a clean slate ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(tidyverse)



#### Read in data ####

geochem.data <- read.csv("dataEdited/geochem/geochem_data_2018.csv",
                         stringsAsFactors = FALSE)


#### Make vector with correct order of samples ####

order.of.samples <- c("2A-P", "2A-N", "2A-K", "2A-A",
                      "3A-O", "3A-K", "3A-F", "3A-A",
                      "WW")
order.of.samples.index <- sapply(order.of.samples,
                                 function(x) {
                                   grep(x, geochem.data$site)
                                 })
geochem.data <- geochem.data[order.of.samples.index, ]


#### Isolate S data for plotting ####

SW_sulfate <- geochem.data %>%
  filter(location == "SW") %>%
  select(site, sulfate_mg_L)
PW_sulfide <- geochem.data %>%
  filter(location == "PW") %>%
  select(site, sulfide_ug_L)



#### Plot PW sulfide and SW sulfate data ####

pdf("results/geochemistry/S_gradients_2018.pdf",
    height = 6,
    width = 8)
par(mfrow = c(2, 1),
    mar = c(3, 3, 3, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
barplot(SW_sulfate$sulfate_mg_L,
        names.arg = SW_sulfate$site,
        ylim = c(0, 60),
        ylab = "Sulfate (mg/L)",
        main = "Surface water sulfate")
barplot(PW_sulfide$sulfide_ug_L,
        names.arg = PW_sulfide$site,
        ylim = c(0, 8000),
        ylab = "Sulfide (µg/L)",
        main = "Porewater sulfide")
dev.off()





#### Isolate Hg data for plotting ####

PW_mehg <- geochem.data %>%
  filter(location == "PW") %>%
  select(site, FMHg_ng_L)
SW_mehg <- geochem.data %>%
  filter(location == "SW") %>%
  select(site, FMHg_ng_L)
SW_thg <- geochem.data %>%
  filter(location == "SW") %>%
  select(site, FTHg_ng_L)
PW_thg <- geochem.data %>%
  filter(location == "PW") %>%
  select(site, FTHg_ng_L)




#### Plot MeHg and Hg at both PW and SW ####

pdf("results/geochemistry/Hg_gradients_2018.pdf",
    width = 12,
    height = 9)
par(mfrow = c(2,2),
    mar = c(3, 3, 3, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
barplot(PW_mehg$FMHg_ng_L,
        names.arg = PW_mehg$site,
        ylim = c(0, 2),
        ylab = "MeHg (ng/L)",
        main = "Porewater MeHg")
barplot(rep(0, length(PW_thg$FTHg_ng_L)),
        names.arg = PW_mehg$site,
        ylim = c(0, 3),
        ylab = "THg (ng/L)",
        main = "No porewater THg data yet")
barplot(SW_mehg$FMHg_ng_L,
        names.arg = PW_mehg$site,
        ylim = c(0, 2),
        ylab = "MeHg (ng/L)",
        main = "Surface water MeHg")
barplot(SW_thg$FTHg_ng_L,
        names.arg = SW_thg$site,
        ylim = c(0, 3.5),
        ylab = "THg (ng/L)",
        main = "Surface water THg")
dev.off()
