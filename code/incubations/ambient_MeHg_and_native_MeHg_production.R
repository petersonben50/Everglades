#### code/incubations/ambient_MeHg_and_native_MeHg_production.R ####
# Benjamin D. Peterson


#### Get cleaned up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(patchwork)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) 


#### Read in porewater data, order it, and generate percent MeHg ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(MeHg_porewater = as.numeric(FMHG),
         HgT_porewater = as.numeric(FTHG),
         percent_MeHg_porewater = MeHg_porewater / HgT_porewater * 100) %>%
  select(siteID, MeHg_porewater, HgT_porewater, percent_MeHg_porewater)


#### Read in hgcA abundance data ####
hgcA.depth <- readRDS("dataEdited/assembly_analysis/hgcA/hgcA_abundance_site.rds") %>%
  group_by(siteID) %>%
  summarise(coverage = mean(coverage))



#### Calculate ambient MeHg in sediment from incubation cores ####
amb.Hg.sediment <- inc.Hg.data %>%
  group_by(coreID) %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  summarise(SMHG_amb_mean = mean(SMHG_amb),
            SMHG_amb_sd = sd(SMHG_amb),
            SMHG_amb_count = n(),
            SMHG_amb_se = SMHG_amb_sd / sqrt(SMHG_amb_count),
            STHG_amb_mean = mean(STHG_amb),
            STHG_amb_sd = sd(STHG_amb),
            STHG_amb_count = n(),
            STHG_amb_se = STHG_amb_sd / sqrt(STHG_amb_count),
            SMHG_amb_percent_mean = mean(SMHG_amb_percent),
            SMHG_amb_percent_sd = sd(SMHG_amb_percent),
            SMHG_amb_percent_count = n(),
            SMHG_amb_percent_se = SMHG_amb_percent_sd / sqrt(SMHG_amb_percent_count)) %>%
  select(coreID, SMHG_amb_mean, SMHG_amb_se,
         STHG_amb_mean, STHG_amb_se,
         SMHG_amb_percent_mean, SMHG_amb_percent_se)
amb.Hg.sediment



#### Plot ambient MeHg in soil ####
MeHg.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_amb_mean - SMHG_amb_se,
                    ymax = SMHG_amb_mean + SMHG_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "MeHg concentrations in sediment",
       y = "MeHg (ng/g)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Linear correlation: ambient MeHg in soil and hgcA ####
amb.Hg.sediment.hgcA <- full_join(amb.Hg.sediment %>% rename(siteID = coreID),
                                  hgcA.depth)
# First check linearity
amb.Hg.sediment.hgcA %>%
  ggplot(aes(x = coverage,
             y = SMHG_amb_mean,
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  theme_bw()
# Actually looks pretty good.

linear.model.ambHg.hgcA <- lm(SMHG_amb_mean ~ coverage,
                              data = amb.Hg.sediment.hgcA)
summary(linear.model.ambHg.hgcA)
# Check residuals
shapiro.test(linear.model.ambHg.hgcA$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.ambHg.hgcA$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.ambHg.hgcA$residuals)
qqline(linear.model.ambHg.hgcA$residuals)
# Not particularly well distributed. Let's try the log transformation
# Although, the limited number of samples might just throw this off. 

# First check linearity
amb.Hg.sediment.hgcA %>%
  ggplot(aes(x = log(coverage, 10),
             y = log(SMHG_amb_mean, 10),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  theme_bw()
# Actually looks pretty good.

linear.model.ambHg.hgcA <- lm(log(SMHG_amb_mean, 10) ~ log(coverage, 10),
                              data = amb.Hg.sediment.hgcA)
# Check residuals
shapiro.test(linear.model.ambHg.hgcA$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.ambHg.hgcA$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.ambHg.hgcA$residuals)
qqline(linear.model.ambHg.hgcA$residuals)
# Large deviation due to a single point (2A-A).
# Overall though, we do have a solid correlation here. 
summary(linear.model.ambHg.hgcA)


#### Plot linear correlation of ambient peat MeHg and hgcA ####

# Calculate p-value
f <- summary(linear.model.ambHg.hgcA)$fstatistic
p.value <- pf(f[1],
              f[2],f[3],
              lower.tail = F) %>%
  round(4)

# Make the plot with log transformed data
ambMeHg.hgcA.plot.logData <- amb.Hg.sediment.hgcA %>%
  ggplot(aes(x = log(coverage, 10),
             y = log(SMHG_amb_mean, 10),
             colour = siteID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  scale_shape_manual(values = point.vector[unique(amb.Hg.sediment.hgcA$siteID)], name = "Site ID") +
  scale_color_manual(values = color.vector[unique(amb.Hg.sediment.hgcA$siteID)], name = "Site ID") +
  scale_fill_manual("black") +
  labs(x = "hgcA abundance (log %)",
       y = "Ambient MeHg in peat (ng/L)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.18, 0.7)) +
  geom_abline(slope = coef(linear.model.ambHg.hgcA)[[2]],
              intercept = coef(linear.model.ambHg.hgcA)[[1]]) +
  geom_label(x = 0.75, y = 1.0,
             label = paste("Adjusted r2 = ", round(summary(linear.model.ambHg.hgcA)$adj.r.squared, 3), "\n",
                           "p = ", p.value,
                           sep = ""),
             color = "black")
ambMeHg.hgcA.plot.logData


#### Plot percent ambient MeHg ####
MeHgPercent.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_percent_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = SMHG_amb_percent_mean - SMHG_amb_percent_se,
                    ymax = SMHG_amb_percent_mean + SMHG_amb_percent_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "Percent MeHg in sediment",
       y = "Percent MeHg") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Plot ambient percent MeHg against hgcA ####
# First check linearity
amb.Hg.sediment.hgcA %>%
  ggplot(aes(x = coverage,
             y = SMHG_amb_percent_mean,
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  theme_bw()
# Check linear model
linear.model.ambHg.hgcA <- lm(SMHG_amb_percent_mean ~ coverage,
                              data = amb.Hg.sediment.hgcA)
summary(linear.model.ambHg.hgcA)


#### Generate plot of HgT levels ####
HgT.soil.plot <- amb.Hg.sediment %>%
  ggplot(aes(x = coreID,
             y = STHG_amb_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = STHG_amb_mean - STHG_amb_se,
                    ymax = STHG_amb_mean + STHG_amb_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(title = "HgT in sediment",
       y = "HgT (ng/g)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Stack plots together ####
MeHg.soil.plot / HgT.soil.plot / MeHgPercent.soil.plot
rm(MeHg.soil.plot, MeHgPercent.soil.plot, HgT.soil.plot)




#### Porewater MeHg concentration plot ####
MeHg.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = MeHg_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "MeHg in porewater",
       y = "MeHg (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater HgT concentration plot ####
HgT.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = HgT_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "HgT in porewater",
       y = "HgT (ng/L)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater MeHg percentage plot ####
MeHg.percent.porewater.plot <- porewater.data %>%
  ggplot(aes(x = siteID,
             y = percent_MeHg_porewater,
             width = 0.8,
             fill = siteID)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = color.vector) +
  labs(title = "Percent MeHg in porewater",
       y = "MeHg (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")


#### Porewater plots stacked ####
MeHg.porewater.plot / HgT.porewater.plot / MeHg.percent.porewater.plot
rm(MeHg.porewater.plot, HgT.porewater.plot, MeHg.percent.porewater.plot)



#### Linear regression: Porewater MeHg percentage and hgcA ####
amb.Hg.porewater.hgcA <- full_join(porewater.data,
                                   hgcA.depth)

# First check linearity
amb.Hg.porewater.hgcA %>%
  ggplot(aes(x = log(coverage),
             y = log(percent_MeHg_porewater),
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  theme_bw()
# Hmm, not particularly linear

linear.model.ambHg.hgcA <- lm(percent_MeHg_porewater ~ coverage,
                              data = amb.Hg.porewater.hgcA)
summary(linear.model.ambHg.hgcA)
# Check residuals
shapiro.test(linear.model.ambHg.hgcA$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.ambHg.hgcA$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.ambHg.hgcA$residuals)
qqline(linear.model.ambHg.hgcA$residuals)
# Not particularly well distributed. Let's try the log transformation
# Although, the limited number of samples might just throw this off. 


#### Calculate: production of MeHg using native porewaters ####
native.porewater.production <- inc.Hg.data %>%
  filter(coreID == matrixID) %>%
  group_by(coreID) %>%
  summarise(meth_spike_per_mean = mean(SMHG_201_percent * 100),
            meth_spike_per_sd = sd(SMHG_201_percent * 100),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count)) %>%
  rename(siteID = coreID) %>%
  select(siteID, meth_spike_per_mean, meth_spike_per_se)
# Add in hgcA data
native.porewater.production.hgcA <- full_join(native.porewater.production,
                                              hgcA.depth)



#### Plot MeHg production using native cores ####
native.porewater.production %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  scale_fill_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       title = "Hg methylation using ambient porewater") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))



#### Native MeHg production vs. MeHg in porewater ####
# Generate data frame for MeHg percent in porewater against MeHg production
MeHg.prod.and.porewater <- native.porewater.production %>%
  left_join(porewater.data %>% rename(coreID = siteID))
# Make plot
MeHg.prod.vs.porewater.MeHg.levels <- MeHg.prod.and.porewater %>%
  ggplot(aes(x = MeHg_porewater,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(x = "MeHg (ng/L)",
       y = "Methylated spike (%)",
       title = "Native Hg methylation vs. porewater MeHg levels") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# Make plot
MeHg.prod.vs.porewater.MeHg.percent <- MeHg.prod.and.porewater %>%
  ggplot(aes(x = percent_MeHg_porewater,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(x = "MeHg (%)",
       y = "Methylated spike (%)",
       title = "Native Hg methylation vs. porewater MeHg percent") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))



#### Native MeHg production vs. MeHg in sediment ####
# MeHg percent in sediment against MeHg production
MeHg.prod.and.sediment <- native.porewater.production %>%
  left_join(amb.Hg.sediment)


# Plot of MeHg production vs. MeHg levels in seds
MeHg.prod.vs.sediment.MeHg.levels <- MeHg.prod.and.sediment %>% 
  ggplot(aes(x = SMHG_amb_mean,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       x = "MeHg levels (ng/g)",
       title = "Native Hg methylation vs. sediment MeHg levels") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# Plot of MeHg production vs. MeHg percent in seds
MeHg.prod.vs.sediment.MeHg.percent <- MeHg.prod.and.sediment %>% 
  ggplot(aes(x = SMHG_amb_percent_mean,
             y = meth_spike_per_mean,
             width = 0.8,
             color = coreID)) + 
  geom_point(aes(size = 3)) +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.001) +
  scale_color_manual(values = color.vector) +
  labs(y = "Methylated spike (%)",
       x = "MeHg (%)",
       title = "Native Hg methylation vs. sediment MeHg percent") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))


#### Plot figures together ####
(MeHg.prod.vs.porewater.MeHg.levels + MeHg.prod.vs.porewater.MeHg.percent) / (MeHg.prod.vs.sediment.MeHg.levels + MeHg.prod.vs.sediment.MeHg.percent)



#### Linear regression: Native MeHg production vs. hgcA ####

# First check linearity
native.porewater.production.hgcA %>%
  ggplot(aes(x = coverage,
             y = meth_spike_per_mean,
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  theme_bw()
# Hmm, not particularly linear

linear.model.ambHg.hgcA <- lm(meth_spike_per_mean ~ coverage,
                              data = native.porewater.production.hgcA)
summary(linear.model.ambHg.hgcA)
# Check residuals
shapiro.test(linear.model.ambHg.hgcA$residuals)
par(mfrow = c(1,2))
plot(density(linear.model.ambHg.hgcA$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(linear.model.ambHg.hgcA$residuals)
qqline(linear.model.ambHg.hgcA$residuals)




#### Run linear regressions of all combos ####

# Native MeHg production vs sediment MeHg concentration
summary(lm(meth_spike_per_mean ~ MeHg_porewater,
           data = MeHg.prod.and.porewater))
# No linear regression

# Native MeHg production vs sediment MeHg percent
summary(lm(meth_spike_per_mean ~ percent_MeHg_porewater,
           data = MeHg.prod.and.porewater))
# No linear regression

# Native MeHg production vs sediment MeHg concentration
summary(lm(meth_spike_per_mean ~ SMHG_amb_mean,
           data = MeHg.prod.and.sediment))
# No linear regression

# Native MeHg production vs sediment MeHg percent
summary(lm(meth_spike_per_mean ~ SMHG_amb_percent_mean,
           data = MeHg.prod.and.sediment))
# No linear regression
