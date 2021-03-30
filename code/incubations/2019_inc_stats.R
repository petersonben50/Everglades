

rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####

MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
PW.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8", "FI", "CYSTEINE", "CONTROL")


#### Prep color and point vector ####
color.vector <- cb.translator[1:6]
names(color.vector) <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
point.vector <- c(9, 18:15, 25)
names(point.vector) <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Read in incubation data and order it ####
rel.methylation.data <- readRDS("dataEdited/2019_incubations/rel_methylation_microbes.rds") %>%
  select(coreID, matrixID, rel_meth_spike) %>%
  rename(siteID = coreID)


#### Read in hgcA abundance data ####
hgcA.depth <- readRDS("dataEdited/2019_analysis_assembly/hgcA/hgcA_abundance_site.rds") %>%
  group_by(siteID) %>%
  summarise(coverage = mean(coverage))
all.data <- rel.methylation.data %>%
  left_join(hgcA.depth)


#### Run a two-way ANOVA on this with hgcA coverage as an independent continuous variable ####
linear.model <- lm(rel_meth_spike ~ coverage,
                   data = all.data)
summary(linear.model)


#### Generate plots ####
all.data %>%
  ggplot(aes(x = coverage,
             y = rel_meth_spike,
             colour = siteID)) +
  geom_point(size = 3,
             aes(shape = siteID)) +
  scale_shape_manual(values = point.vector) +
  scale_color_manual(values = color.vector) +
  scale_fill_manual("black") +
  labs(x = "Read coverage of hgcA sequences",
       y = "Relative methylation potential\nof sediment cores",
       title = "Relative methylation potential vs.hgcA coverage in sediment") +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black")) +
  geom_smooth(method = "lm",
              formula = y ~ x)
