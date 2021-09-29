#### code/incubations/RMP_analysis_porewater.R ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
source("code/setup_PW_core_order_color_points.R")


#### Read in porewater data ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis for Brett and Ben.xlsx",
                            sheet = "porewater_data") %>%
  mutate(siteID = fct_relevel(siteID, MG.order),
         DOC = gsub("\\*", "", DOC_mg.L) %>% as.numeric(),
         SUVA = gsub("\\*", "", `SUVA_L.mg*m`) %>% as.numeric(),
         UV_254_AU = as.numeric(`UV_254_AU.cm-1`),
         sulfate_mg.L = gsub("nd", 0, sulfate_mg.L) %>% as.numeric(),
         sulfide_µg.L = gsub("nd", 0, sulfide_µg.L) %>% as.numeric()) %>%
  rename(matrixID = siteID) %>%
  select(matrixID, DOC, SUVA, UV_254_AU, sulfate_mg.L, sulfide_µg.L)


#### Read in RMP data ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Join RMP and porewater data ####
Hg.inc.data <- Hg.inc.data %>%
  select(coreID, matrixID, RMP_porewater) %>%
  rename(siteID = coreID)
all.data <- Hg.inc.data %>%
  left_join(porewater.data)


#### RMP of porewaters along sulfate gradient ####
all.data %>%
  group_by(matrixID, siteID) %>%
  summarise(RMP_porewater = mean(RMP_porewater)) %>%
  ggplot(aes(x = matrixID,
             y = RMP_porewater,
             group = siteID,
             color = siteID)) +
  geom_point() +
  geom_line() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = c(0.9, 0.75))



#### RMP of porewaters against environmental variables ####
RMP.sulfide <- all.data %>%
  ggplot(aes(x = sulfide_µg.L,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Sulfide (µg/L)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = c(0.6, 0.35))
RMP.sulfide
# Sulfide spans too many orders of magnitude to accurately visualize this.
# Let's do a log transformation of the sulfide
RMP.sulfide.log <- all.data %>%
  ggplot(aes(x = log(sulfide_µg.L, 10),
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Log sulfide (µg/L)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = c(0.7, 0.45))

RMP.SUVA <- all.data %>%
  filter(!is.na(SUVA)) %>%
  ggplot(aes(x = SUVA,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "SUVA (L/mg/m)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = "none")

RMP.DOC <- all.data %>%
  filter(!is.na(DOC)) %>%
  ggplot(aes(x = DOC,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "DOC (mg/L)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = "none")


RMP.UV <- all.data %>%
  mutate(UV_254_AU = DOC * SUVA / 100) %>%
  filter(!is.na(UV_254_AU)) %>%
  ggplot(aes(x = UV_254_AU,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Site ID",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "UV 254 (cm-1)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = "none")



#### Plot em all out ####
(RMP.sulfide.log + RMP.SUVA) / (RMP.DOC + RMP.UV)



#### Prettify the RMP vs. sulfide plot ####
RMP.sulfide.log.for.saving <- all.data %>%
  ggplot(aes(x = log(sulfide_µg.L, 10),
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Log sulfide (µg/L)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.7, 0.45))
# As RDS
saveRDS(object = RMP.sulfide.log.for.saving,
        file = "results/incubations/RMP_porewater_sulfide.rds")




#### Generate and check linear model of RMP porewater vs. SUVA ####
PW.linear.model <- lm(RMP_porewater ~ SUVA,
                   data = all.data)
# Check residuals
par(mfrow = c(1,2))
plot(density(PW.linear.model$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(PW.linear.model$residuals)
qqline(PW.linear.model$residuals)
# These have a slight left skew

shapiro.test(PW.linear.model$residuals)
# Normally distributed, statistically speaking!

summary(PW.linear.model)
# Major effect of SUVA on RMP
# Adjusted R-squared: 0.4939
# p-value: 5.593e-14


#### Generate plot with linear model plotted on it ####
rSquared <- round(summary(PW.linear.model)$adj.r.squared, 2)
RMP.SUVA.plot <- all.data %>%
  filter(!is.na(SUVA)) %>%
  ggplot(aes(x = SUVA,
             y = RMP_porewater,
             colour = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  scale_fill_manual("black") +
  labs(x = "SUVA (L/mg/m)",
       y = "Porewater RMP (%)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.82, 0.33)) +
  geom_abline(slope = coef(PW.linear.model)[[2]],
              intercept = coef(PW.linear.model)[[1]]) +
  geom_label(x = 4, y = 82,
             label = paste("Adjusted r2 = ", round(summary(PW.linear.model)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
RMP.SUVA.plot


#### Save out plot of porewater RMP vs SUVA ####
pdf("results/incubations/RMP_porewater_SUVA.pdf",
    height = 5.5,
    width = 7)
RMP.SUVA.plot
dev.off()
# As RDS
saveRDS(object = RMP.SUVA.plot,
        file = "results/incubations/RMP_porewater_SUVA.rds")