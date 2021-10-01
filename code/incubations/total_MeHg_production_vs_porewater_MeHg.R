#### code/incubations ####
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
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")



#### Read in incubation data ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Read in porewater data, order it, and generate percent MeHg ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  filter(siteID %in% MG.order) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(MeHg_porewater = as.numeric(FMHG),
         HgT_porewater = as.numeric(FTHG),
         percent_MeHg_porewater = MeHg_porewater / HgT_porewater * 100) %>%
  select(siteID, MeHg_porewater, HgT_porewater, percent_MeHg_porewater)


#### Calculate mean MeHg production at a site, regardless of porewater used ####
mean.MeHg.production.data <- Hg.inc.data %>%
  group_by(coreID) %>%
  mutate(SMHG_201_percent = SMHG_201_percent * 100) %>%
  summarise(SMHG_201_mean = mean(SMHG_201_percent),
            SMHG_201_sd = sd(SMHG_201_percent),
            SMHG_201_count = n(),
            SMHG_201_se = SMHG_201_sd / sqrt(SMHG_201_count)) %>%
  ungroup() %>%
  rename(siteID = coreID) %>%
  select(siteID, SMHG_201_mean, SMHG_201_sd, SMHG_201_se)


#### Combine data ####
all.data.for.curve <- porewater.data %>%
  full_join(mean.MeHg.production.data)


#### Plot data ####
pdf("results/incubations/total_MeHg_production_vs_MeHg_porewater.pdf",
    width = 6,
    height = 6)
par(mar = c(4.5, 4.5, 1, 1))
plot(y = all.data.for.curve$SMHG_201_mean,
     x = all.data.for.curve$percent_MeHg_porewater,
     col = color.vector[all.data.for.curve$siteID],
     pch = 16,
     xlab = "Percent MeHg in porewater",
     ylab = "MeHg production in core",
     xlim = c(0, 14),
     ylim = c(0, 3))
a <- all.data.for.curve
arrows(a$percent_MeHg_porewater, a$SMHG_201_mean-a$SMHG_201_se,
       a$percent_MeHg_porewater, a$SMHG_201_mean+a$SMHG_201_se,
       code=3,
       length=0.02,
       angle = 90)
legend("bottomright",
       legend = all.data.for.curve$siteID,
       pch = 16,
       col = color.vector[all.data.for.curve$siteID])
#### Calculate log-fitting line ####
model <- lm(SMHG_201_mean ~ log(percent_MeHg_porewater),
            data = all.data.for.curve)
summary(model)
# It is significant!


#### Plot model with confidence band ####
matlines(x=seq(from=1,to=20,length.out=1000),
         y=predict(model, newdata=list(percent_MeHg_porewater=seq(from=1,to=20,length.out=1000)),
                   interval="confidence"))
dev.off()
