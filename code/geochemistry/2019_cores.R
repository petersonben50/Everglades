rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)




#### Read in incubation data ####

inc.Hg.data <- read_xlsx("dataEdited/geochem/incubation_cores_2019_from_Dave_edited.xlsx",
                     sheet = "inc_data_cleaned") %>%
  mutate(siteID = gsub(" WI TRANSECT 1-", "", siteID),
         siteID = gsub("LOX8", "LOX", siteID),
         PW_origin = gsub(" PW", "", PW_origin)) %>%
  filter(siteID == PW_origin) %>%
  select(-PW_origin)


Hg.data <- read_xlsx("dataEdited/geochem/incubation_cores_2019_from_Dave_edited.xlsx",
                     sheet = "Hg_data") %>%
  mutate(siteID = gsub(" WI TRANSECT 1-", "", siteID),
         siteID = gsub(" WI TRANSECT 2-", "", siteID),
         siteID = gsub("LOX8", "LOX", siteID))




#### Order of samples to be plotted ####

sample.order.vector <- 1:6
names(sample.order.vector) <- c("2AN", "2AA", "3AO", "3AN", "3AF", "LOX")

sample.order.index <- sapply(names(sample.order.vector),
                             function(x) {
                               grep(x, Hg.data$siteID)
                               })
Hg.data <- Hg.data[sample.order.index, ]
rm(sample.order.index)

#### Generate vector to order sites ####



#### Compare incubation values to porewater MeHg ####

par(mfrow = c(3, 2),
    mar = c(3, 3, 3, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)

# Add porewater MeHg data
barplot(Hg.data$PW_FMHG,
        names.arg = Hg.data$siteID,
        ylim = c(0, 0.1),
        ylab = "PW MeHg",
        main = "Porewater MeHg")
barplot(Hg.data$PW_FMHG / Hg.data$PW_FTHG,
        names.arg = Hg.data$siteID,
        ylim = c(0, 0.12),
        ylab = "PW MeHg (%)",
        main = "Porewater Percent MeHg")

# Add ambient MeHg levels from cores
plot(x = sample.order.vector[inc.Hg.data$siteID],
     y = inc.Hg.data$SMHG_amb,
     xlab = "Site ID",
     ylab = "MeHg",
     pch = 19,
     xaxt = "n",
     main = "Ambient MeHg in cores")
axis(1,
     at = 1:6,
     labels = names(sample.order.vector))

plot(x = sample.order.vector[inc.Hg.data$siteID],
     y = inc.Hg.data$SMHG_amb / inc.Hg.data$STHG_amb * 100,
     xlab = "Site ID",
     ylab = "MeHg (%)",
     pch = 19,
     xaxt = "n",
     main = "Ambient Percent MeHg in cores")
axis(1,
     at = 1:6,
     labels = names(sample.order.vector))

# Add Me201Hg production data
plot(x = sample.order.vector[inc.Hg.data$siteID],
     y = inc.Hg.data$SMHG_201,
     xlab = "Site ID",
     ylab = "Me201Hg",
     pch = 19,
     xaxt = "n",
     main = "Excess Me201Hg in cores")
axis(1,
     at = 1:6,
     labels = names(sample.order.vector))

plot(x = sample.order.vector[inc.Hg.data$siteID],
     y = inc.Hg.data$SMHG_201 / inc.Hg.data$STHG_201 * 100,
     xlab = "Site ID",
     ylab = "Me201Hg (%)",
     pch = 19,
     xaxt = "n",
     main = "Percent excess Me201Hg in cores")
axis(1,
     at = 1:6,
     labels = names(sample.order.vector))




#### Plot MeHg produced against ambient MeHg in cores ####

par(mfrow = c(1, 2))
plot(x = inc.Hg.data$SMHG_amb,
     y = inc.Hg.data$SMHG_201,
     xlab = "Ambient MeHg",
     ylab = "Me201Hg",
     cex = 0)
text(x = inc.Hg.data$SMHG_amb,
     y = inc.Hg.data$SMHG_201,
     labels = inc.Hg.data$siteID)

plot(x = inc.Hg.data$SMHG_amb / inc.Hg.data$STHG_amb * 100,
     y = inc.Hg.data$SMHG_201 / inc.Hg.data$STHG_201 * 100,
     xlab = "Ambient % MeHg",
     ylab = "% Me201Hg",
     cex = 0)
text(x = inc.Hg.data$SMHG_amb / inc.Hg.data$STHG_amb * 100,
     y = inc.Hg.data$SMHG_201 / inc.Hg.data$STHG_201 * 100,
     labels = inc.Hg.data$siteID)
