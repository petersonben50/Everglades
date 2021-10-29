#### code/incubations/statistical_analysis_microbialCommunity.R ####
# Benjamin D. Peterson

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in hgcA abundance data ####
hgcA.depth <- readRDS("dataEdited/assembly_analysis/hgcA/hgcA_abundance_site.rds") %>%
  group_by(siteID) %>%
  summarise(coverage = mean(coverage))


#### Read in methylation data ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Join methylation and hgcA data ####
Hg.inc.data <- Hg.inc.data %>%
  select(coreID, matrixID, RMP_core) %>%
  rename(siteID = coreID)
all.data <- Hg.inc.data %>%
  left_join(hgcA.depth)


#### Run a linear model on this with hgcA coverage as an independent continuous variable ####
# Check with non-log-transformed data
shapiro.test(all.data$RMP_core)
par(mfrow = c(1,1))
hist(all.data$RMP_core)
# Not normal, let's check residuals
linear.model <- lm(RMP_core ~ coverage,
                   data = all.data)
# Check residuals
hist(linear.model$residuals)
shapiro.test(linear.model$residuals)
# Not bad, but not particularly good either

# Try with log-transformed data
shapiro.test(log(all.data$RMP_core))
hist(log(all.data$RMP_core))
# Still not normal, but doesn't necessarily need to be for lm.
# Generate linear model
linear.model.log <- lm(log(RMP_core) ~ coverage,
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
# Slight left skew

# Do we see linearity?
all.data %>%
  ggplot(aes(x = coverage,
             y = log(RMP_core),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID))
# Not really


#### Linear model with a log-log transformation ####
# Let's check linearity first this time
all.data %>%
  ggplot(aes(x = log(coverage, 10),
             y = log(RMP_core, 10),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID))
# Better. Try the linear regression on
# the log-log transformation
linear.model.log.log <- lm(log(RMP_core, 10) ~ log(coverage, 10),
                           data = all.data)
# Check residuals
shapiro.test(linear.model.log.log$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.log.log$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.log.log$residuals)
qqline(linear.model.log.log$residuals)
# Little bit further from normality, but fairly close.
# Mostly has a right skew here.
# We'll go ahead with it.

# Summarize model
summary(linear.model.log.log)
# Linear fit is even better, 0.4935



#### Generate plot of RMP vs. hgcA ####

# Linear model to use
linear.model.to.use <- linear.model.log.log

# Calculate p-value
f <- summary(linear.model.to.use)$fstatistic
p.value <- pf(f[1],
              f[2],f[3],
              lower.tail = F) %>%
  round(4)

# Make the plot with log transformed data
RMP.hgcA.plot.logData <- all.data %>%
  ggplot(aes(x = log(coverage, 10),
             y = log(RMP_core, 10),
             colour = siteID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  scale_shape_manual(values = point.vector[unique(all.data$siteID)], name = "Site ID") +
  scale_color_manual(values = color.vector[unique(all.data$siteID)], name = "Site ID") +
  scale_fill_manual("black") +
  labs(x = "hgcA abundance (log %)",
       y = "Sediment cores RMP (log %)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.18, 0.7)) +
  geom_abline(slope = coef(linear.model.to.use)[[2]],
              intercept = coef(linear.model.to.use)[[1]]) +
  geom_label(x = 0.75, y = 1.7,
             label = paste("Adjusted r2 = ", round(summary(linear.model.to.use)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")


#### Generate plot of RMP vs. hgcA ####

# Make the plot with log transformed scales, not data
RMP.hgcA.plot <- all.data %>%
  ggplot(aes(x = coverage,
             y = RMP_core,
             colour = siteID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  scale_shape_manual(values = point.vector[unique(all.data$siteID)], name = "Site ID") +
  scale_color_manual(values = color.vector[unique(all.data$siteID)], name = "Site ID") +
  scale_fill_manual("black") +
  scale_x_continuous(trans = 'log10',
                     limits = c(1, 15)) +
  scale_y_continuous(trans = 'log10',
                     limits = c(1, 100)) +
  labs(x = "hgcA abundance (%)",
       y = "Sediment cores RMP (%)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.18, 0.7)) +
  # geom_abline(slope = coef(linear.model.to.use)[[2]],
  #             intercept = coef(linear.model.to.use)[[1]]) +
  geom_label(x = 0.75, y = 1.7,
             label = paste("Adjusted r2 = ", round(summary(linear.model.to.use)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
RMP.hgcA.plot



#### Save out plot ####
pdf("results/incubations/RMP_microbialCommunity_hgcA.pdf",
    height = 5.5,
    width = 7)
RMP.hgcA.plot
dev.off()
saveRDS(object = RMP.hgcA.plot,
        file = "results/incubations/RMP_microbes_hgcA.rds")
