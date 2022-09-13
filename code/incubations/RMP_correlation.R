#### code/incubations/RMP_correlation.R ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in data ####
all.data <- readRDS("dataEdited/incubations/RMP_ambient_MeHg_production_comparison.rds")


#### Plot RMP pore waters vs. RMP core ####
all.data %>%
  ggplot(aes(x = RMP_core_mean,
             y = RMP_porewater_mean,
             col = siteID,
             shape = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector) +
  scale_shape_manual(values = point.vector) +
  theme_bw()


#### Linear regression: RMPcore vs RMPporewater ####
# Try with log-transformed data
# Generate linear model
linear.model.log <- lm(log(RMP_core_mean) ~ log(RMP_porewater_mean),
                       data = all.data)
# Check residuals
shapiro.test(linear.model.log$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.log$residuals)
qqline(linear.model.log$residuals)
# These are not quite normally distributed, but a little bit better
summary(linear.model.log)
