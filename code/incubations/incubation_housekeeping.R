#### code/incubations/incubation_housekeeping.R ####
# Benjamin D. Peterson

# Additional notes in the corresponding Obsidian document:
# Statistical analysis of MeHg production - re-analyzed for paper

# Results here: results/incubation_housekeeping.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggpubr)
library(lme4)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")



#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                         sheet = "cleaned_incubation_data",) %>%
  mutate(coreID = fct_relevel(coreID, MG.order),
         matrixID = fct_relevel(matrixID, PW.order),
         SMHG_amb_percent = SMHG_amb_percent * 100,
         SMHG_201_percent = SMHG_201_percent * 100)



#### Plot: soil variables by core origination ####
LOI.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = LOI,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector,
                     name = "Porewater ID") +
  ylab("Loss on ignition (%)") +
  ylim(c(50, 100)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = c(0.82, 0.33),
        axis.title.x = element_blank())
dryWt.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = dryWt,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector) +
  ylab("Dry weight (%)") +
  ylim(c(0, 30)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none",
        axis.title.x = element_blank())

pdf("results/incubations/LOI_dry_weight.pdf",
    height = 5,
    width = 10)
ggarrange(dryWt.plot,
          LOI.plot,
          labels = c("A.", "B."))
dev.off()


#### Stats: dry weight ####
# Dry weight stats
max(inc.Hg.data$dryWt)
min(inc.Hg.data$dryWt)

inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(dryWt_mean = mean(dryWt),
            dryWt_sd = sd(dryWt),
            dryWt_rsd = dryWt_sd / dryWt_mean * 100,
            max_dryWt = max(dryWt),
            min_dryWt = min(dryWt),
            range_dryWt = max_dryWt - min_dryWt)

#### Stats: LOI ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(LOI_mean = mean(LOI),
            LOI_sd = sd(LOI),
            LOI_rsd = LOI_sd / LOI_mean * 100,
            max_LOI = max(LOI),
            min_LOI = min(LOI),
            range_LOI = max_LOI - min_LOI)



#### Stats: Check relationship between these soil parameters and sediment core/porewater source ####
# This is mostly paranoia.

# LOI
twoWayAnova_loi <- aov(LOI ~ coreID * matrixID,
                       data = inc.Hg.data)
par(mfrow = c(1,2))
hist(inc.Hg.data$LOI)
hist(log(inc.Hg.data$LOI))
# Check residuals
shapiro.test(twoWayAnova_loi$residuals)
qqnorm(twoWayAnova_loi$residuals)
qqline(twoWayAnova_loi$residuals)
plot(density(twoWayAnova_loi$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# Residuals are fucked with or without the log transformation.
# Let's remove it.
summary(twoWayAnova_loi)
# No effect of porewater source on LOI, but there is an effect of the core, of course
# No interaction either

# Percent dry weight
twoWayAnova_dryWt <- aov(dryWt ~ coreID * matrixID,
                         data = inc.Hg.data)
par(mfrow = c(1,2))
hist(inc.Hg.data$dryWt)
hist(log(inc.Hg.data$dryWt))
shapiro.test(twoWayAnova_dryWt$residuals)
# Normal distribution of residuals
qqnorm(twoWayAnova_dryWt$residuals)
qqline(twoWayAnova_dryWt$residuals)
plot(density(twoWayAnova_dryWt$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# Not perfect, but somewhat okay.
# This is actually better without the log transformation, so let's remove that.
summary(twoWayAnova_dryWt)
# Both porewater source and core source have effect on dry weight...
# Now why would the porewater source influence that? Not sure that it's
# worth worrying about, but something to keep in the back of my mind.


#### Plot: Porewater source effects on dry weight ####
inc.Hg.data %>%
  ggplot(aes(x = matrixID,
             y = dryWt)) +
  geom_boxplot() +
  geom_point(aes(col = coreID),
             position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector) +
  ylab("Dry weight (%)") +
  theme_bw()



#### Plot: Ambient Hg variables by site ####
# MeHg plot
MeHg.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector,
                     name = "Porewater origin") +
  ylim(c(0, 9)) +
  ylab("MeHg (ng/g)") +
  xlab("Peat core source") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")

# HgT plot
HgT.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = STHG_amb,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector) +
  ylim(c(0, 400)) +
  ylab("HgT (ng/g)") +
  xlab("Peat core source") +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        legend.position = "none")

# Percent MeHg plot
MeHgPercent.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = SMHG_amb_percent,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector) +
  ylab("Percent MeHg") +
  theme_bw() +
  theme(legend.position = "none")

MeHg.plot + HgT.plot + MeHgPercent.plot


#### Plot: Amended Hg isotope concentrations by site ####

# Concentration of 201HgT plot
iso.HgT.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = STHG_201,
             col = matrixID)) +
  geom_point(position = position_dodge(width=0.3)) +
  scale_color_manual(values = color.vector) +
  ylab(bquote(""^201~HgT ~~(ng/g))) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(c(0, 70))
# Concentration of Me201Hg plot
iso.MeHg.plot <- inc.Hg.data %>%
  ggplot(aes(x = coreID,
             y = SMHG_201,
             col = matrixID)) +
  geom_jitter(width = 0.1) +
  scale_color_manual(values = color.vector) +
  ylab(bquote(~Me^201~Hg ~~(ng/g))) +
  theme_bw() +
  theme(legend.position = "none") +
  ylim(c(0, 3.8))


#### Saved plot: All measured ambient constituents ####
pdf("results/incubations/ambient_measurements.pdf",
    height = 9,
    width = 7.2)
ggarrange(dryWt.plot, LOI.plot,
          HgT.plot, MeHg.plot,
          iso.HgT.plot, iso.MeHg.plot,
          nrow = 3, ncol = 2,
          labels = c("a", "b", "c", "d", "e", "f"))
dev.off()





#### Stats: HgT ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(HGT_mean = mean(STHG_amb),
            HGT_sd = sd(STHG_amb),
            HGT_rsd = HGT_sd / HGT_mean * 100,
            max_HGT = max(STHG_amb),
            min_HGT = min(STHG_amb),
            range_HGT = max_HGT - min_HGT)


#### Stats: MeHg ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(MEHG_mean = mean(SMHG_amb),
            MEHG_sd = sd(SMHG_amb),
            MEHG_rsd = MEHG_sd / MEHG_mean * 100,
            max_MEHG = max(SMHG_amb),
            min_MEHG = min(SMHG_amb),
            range_MEHG = max_MEHG - min_MEHG)



#### Stats: Amended 201HgT ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(HGT_mean = mean(STHG_201),
            HGT_sd = sd(STHG_201),
            HGT_rsd = HGT_sd / HGT_mean * 100,
            max_HGT = max(STHG_201),
            min_HGT = min(STHG_201),
            range_HGT = max_HGT - min_HGT)



#### Normalized ambient MeHg data by core ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(SMHG_amb_percent_max = max(SMHG_amb_percent)) %>%
  right_join(inc.Hg.data %>% select(coreID, matrixID, SMHG_amb_percent)) %>%
  left_join(inc.Hg.data %>%
              group_by(coreID, matrixID) %>%
              summarise(site_PW_mean = mean(SMHG_amb_percent))) %>%
  mutate(percent_shift_MeHg = SMHG_amb_percent / SMHG_amb_percent_max * 100) %>%
  ggplot(aes(x = coreID,
             y = percent_shift_MeHg,
             color = matrixID)) +
  geom_jitter(width = 0.1) +
  ylab("Percent shift in MeHg relative to maximum in that group") +
  scale_color_manual(values = color.vector) +
  theme_bw()



#### Normalized ambient MeHg data by porewater ####
inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(SMHG_amb_percent_max = max(SMHG_amb_percent)) %>%
  right_join(inc.Hg.data %>% select(coreID, matrixID, SMHG_amb_percent)) %>%
  left_join(inc.Hg.data %>%
              group_by(coreID) %>%
              summarise(SMHG_amb_percent_max = max(SMHG_amb_percent)) %>%
              right_join(inc.Hg.data %>% select(coreID, matrixID, SMHG_amb_percent)) %>%
              mutate(percent_shift_MeHg = SMHG_amb_percent / SMHG_amb_percent_max * 100) %>%
              group_by(coreID, matrixID) %>%
              summarise(site_mean_percent_shift_MeHg_percent = mean(percent_shift_MeHg))) %>%
  mutate(percent_shift_MeHg = SMHG_amb_percent / SMHG_amb_percent_max * 100) %>%
  ggplot(aes(x = matrixID,
             y = percent_shift_MeHg)) +
  geom_boxplot() +
  geom_line(aes(y = site_mean_percent_shift_MeHg_percent,
                group = coreID,
                color = coreID)) +
  # geom_jitter(aes(color = coreID),
  #             width = 0.1) +
  ylab("Normalized shift in MeHg percent\nrelative to maximum in that group") +
  scale_color_manual(values = color.vector) +
  theme_bw()




#### Check relationship between MeHg constituents and sediment core and porewater source ####

#### Effects of porewater/sediment on MeHg levels ####
# Do log transform, the Q-Q plot was pretty bad
# when I didn't do the log-transform
twoWayAnova_MeHg <- aov(log(SMHG_amb, 10) ~ coreID * matrixID,
                       data = inc.Hg.data)
# Check out residuals
par(mar = c(4.5, 4.5, 3, 1),
    mfrow = c(1,2))
shapiro.test(twoWayAnova_MeHg$residuals)
qqnorm(twoWayAnova_MeHg$residuals)
qqline(twoWayAnova_MeHg$residuals)
plot(density(twoWayAnova_MeHg$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# Distribution of residuals is pretty good with the log transformation.
summary(twoWayAnova_MeHg)
# There is barely an effect of the porewater on the ambient
# MeHg levels (p = 0.049). Especially tiny compared to effect
# of sediment core on MeHg.
# No interactive effect.



#### Effects of porewater/sediment on HgT levels ####
# Marginal improvement of residual distribution when
# data is log transformed
twoWayAnova_HgT <- aov(log(STHG_amb) ~ coreID * matrixID,
                         data = inc.Hg.data)
par(mar = c(4.5, 4.5, 3, 1),
    mfrow = c(1,2))
shapiro.test(twoWayAnova_HgT$residuals)
qqnorm(twoWayAnova_HgT$residuals)
qqline(twoWayAnova_HgT$residuals)
plot(density(twoWayAnova_HgT$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# Pretty good distribution of residuals again.
summary(twoWayAnova_HgT)
# Huge effect of core source, but no effect of porewater origin.


#### Effects of porewater/sediment on percent MeHg ####
twoWayAnova_MeHgPercent_raw <- aov(SMHG_amb_percent ~ coreID * matrixID,
                               data = inc.Hg.data)
twoWayAnova_MeHgPercent <- aov(log(SMHG_amb_percent) ~ coreID * matrixID,
                               data = inc.Hg.data)
hist(twoWayAnova_MeHgPercent_raw$residuals,
     breaks = 20)
hist(twoWayAnova_MeHgPercent$residuals,
     breaks = 20)
shapiro.test(twoWayAnova_MeHgPercent$residuals)
# Need log transform to get residuals to normal
qqnorm(twoWayAnova_MeHgPercent$residuals)
qqline(twoWayAnova_MeHgPercent$residuals)
plot(density(twoWayAnova_MeHgPercent$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
summary(twoWayAnova_MeHgPercent)
# Major effect of core ID. No effect of porewater matrix.
# Very very small interactive effect.
