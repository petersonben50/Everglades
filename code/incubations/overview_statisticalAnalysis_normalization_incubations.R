#### code/incubations/overview_statisticalAnalysis_normalization_incubations.R ####
# Benjamin D. Peterson

# Additional notes in the corresponding Obsidian document:
# Overview, statistical analysis, normalization of incubations

# Powerpoint with notes: results/overview_statisticalAnalysis_normalization_incubations

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order),
         matrixID = fct_relevel(matrixID, PW.order))


#### 




#### Test for normality ####
# Normality of overall data set
shapiro.test(inc.Hg.data$SMHG_201_percent)
hist(inc.Hg.data$SMHG_201_percent,
     breaks = 20)
# Not even close to normal.

# Log transformed data
shapiro.test(log(inc.Hg.data$SMHG_201_percent, 10))
hist(log(inc.Hg.data$SMHG_201_percent, 10),
     breaks = 20)
# Still not close.


#### Check validity of ANOVA as analytical tool ####

# Let's see if the residuals of an ANOVA are normally
# distributed, since that seems to be the main issue.

# Non-transformed data
twoWayAnova_inc <- aov(SMHG_201_percent ~ coreID * matrixID,
                       data = inc.Hg.data)
shapiro.test(twoWayAnova_inc$residuals)
hist(twoWayAnova_inc$residuals,
     breaks = 40)
qqnorm(twoWayAnova_inc$residuals)
# Looks bad...
rm(twoWayAnova_inc)


# Log-transformed data
twoWayAnova_inc_log <- aov(log(SMHG_201_percent, 10) ~ coreID * matrixID,
                           data = inc.Hg.data)
par(mfrow = c(1,2))
hist(twoWayAnova_inc_log$residuals,
     breaks = 20)
plot(density(twoWayAnova_inc_log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
shapiro.test(twoWayAnova_inc_log$residuals)
# QQ-normal plot
qqnorm(twoWayAnova_inc_log$residuals)
qqline(twoWayAnova_inc_log$residuals)




# Square-root-transformed data
twoWayAnova_inc_sqrt <- aov(sqrt(SMHG_201_percent) ~ coreID * matrixID,
                            data = inc.Hg.data)
par(mfrow = c(1,2))
hist(twoWayAnova_inc_sqrt$residuals,
     breaks = 20)
plot(density(twoWayAnova_inc_sqrt$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
shapiro.test(twoWayAnova_inc_sqrt$residuals)
# QQ-normal plot
qqnorm(twoWayAnova_inc_sqrt$residuals)
qqline(twoWayAnova_inc_sqrt$residuals)




#### Storming ahead with two-way ANOVA ####
summary(twoWayAnova_inc_log)

