#### code/incubations/2019_inc_stats.R ####
# Benjamin D. Peterson

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
# PW.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8", "FI", "CYSTEINE", "CONTROL")


#### Prep color and point vector ####
color.vector <- cb.translator[1:6]
names(color.vector) <- MG.order
point.vector <- c(9, 18:15, 25)
names(point.vector) <- MG.order



#### Read in incubation data and order it ####
rel.methylation.data <- readRDS("dataEdited/2019_incubations/rel_methylation_microbes.rds") %>%
  select(coreID, matrixID, MeHg_fract_spike, rel_meth_spike) %>%
  filter(matrixID %in% names(color.vector))



#### Test whether or not both the PW and the core have significant effects on MeHg production ####
# Normality of raw data
shapiro.test(rel.methylation.data$MeHg_fract_spike)
hist(rel.methylation.data$MeHg_fract_spike)
# Not even close to normal.

# Log transformed data
shapiro.test(rel.methylation.data$MeHg_fract_spike)
hist(log(rel.methylation.data$MeHg_fract_spike, 10))
# Still not really normal... 

# Let's see if the residuals of an ANOVA are normally
# distributed, since that seems to be the main issue.
# We'll do it with both
twoWayAnova_inc <- aov(log(MeHg_fract_spike, 10) ~ coreID * matrixID,
                       data = rel.methylation.data)
# twoWayAnova_inc <- aov(MeHg_fract_spike ~ coreID * matrixID,
#                        data = rel.methylation.data)
shapiro.test(twoWayAnova_inc$residuals)
hist(twoWayAnova_inc$residuals,breaks = 20)
# The log-transformed data gets close, but not quite there.
# Let's check it anyways:
summary(twoWayAnova_inc)


#### Read in hgcA abundance data ####
hgcA.depth <- readRDS("dataEdited/2019_analysis_assembly/hgcA/hgcA_abundance_site.rds") %>%
  group_by(siteID) %>%
  summarise(coverage = mean(coverage))


#### Join methylation and hgcA data ####
rel.methylation.data <- rel.methylation.data %>%
  select(coreID, matrixID, rel_meth_spike) %>%
  rename(siteID = coreID)
all.data <- rel.methylation.data %>%
  left_join(hgcA.depth)


#### Run a linear model on this with hgcA coverage as an independent continuous variable ####

# Check with non-log-transformed data
shapiro.test(rel.methylation.data$rel_meth_spike)
hist(rel.methylation.data$rel_meth_spike)
# Not normal, let's check residuals
linear.model <- lm(rel_meth_spike ~ coverage,
                   data = all.data)
summary(linear.model)
# Check residuals
hist(linear.model$residuals)
shapiro.test(linear.model$residuals)

# Try with log-transformed data
shapiro.test(log(rel.methylation.data$rel_meth_spike))
hist(log(rel.methylation.data$rel_meth_spike))
# Still not normal, but doesn't necessarily need to be for lm.
# Generate linear model
linear.model.log <- lm(log(rel_meth_spike) ~ coverage,
                       data = all.data)
summary(linear.model.log)
# Check residuals
hist(linear.model.log$residuals)
shapiro.test(linear.model.log$residuals)
# These are normally distributed
# Do we see linearity?
all.data %>%
  ggplot(aes(x = coverage,
             y = log(rel_meth_spike),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID))
# Not really. 

# Let's check out a  log 
# transformation of both
# variables.
all.data %>%
  ggplot(aes(x = log(coverage),
             y = log(rel_meth_spike),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID))
# Better. Try the linear regression on
# the log-log transformation
linear.model.log.log <- lm(log(rel_meth_spike) ~ log(coverage),
                           data = all.data)
summary(linear.model.log.log)
# Linear fit is even better, 0.527
# Check residuals
hist(linear.model.log.log$residuals)
shapiro.test(linear.model.log.log$residuals)
# Lost normality of the residuals there. Hmm. 
# I think we still need to go with this. 


#### Generate plots ####

# Linear model to use
linear.model.to.use <- linear.model.log.log

# Calculate p-value
f <- summary(linear.model.to.use)$fstatistic
p.value <- pf(f[1],
              f[2],f[3],
              lower.tail = F) %>%
  round(4)

# Make the plot
log.log.plot <- all.data %>%
  ggplot(aes(x = log(coverage),
             y = log(rel_meth_spike),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  scale_shape_manual(values = point.vector) +
  scale_color_manual(values = color.vector) +
  scale_fill_manual("black") +
  labs(x = "Log abundance of hgcA ",
       y = "Log RMP of sediment cores",
       title = element_blank()) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = "none") +
  geom_abline(slope = coef(linear.model.to.use)[[2]],
              intercept = coef(linear.model.to.use)[[1]]) +
  geom_label(x = 2.2, y = -3.2,
             label = paste("Adjusted r2 = ", round(summary(linear.model.to.use)$adj.r.squared, 2), "\n",
                           "p < 0.0001",
                           sep = ""),
             color = "black")


#### Save out PDF ####
pdf("results/2019_incubations/MeHg_vs_hgcA_LM.pdf",
    width = 5,
    height = 4)
log.log.plot
dev.off()
