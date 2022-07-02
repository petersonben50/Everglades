#### code/incubations/RMP_analysis_porewater.R ####
# Benjamin D. Peterson

# Obsidian notes here: Incubations - statistical analysis porewaters
# Results: results/RMP_porewater.pptx

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lme4)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Data loading: Porewater data ####
porewater.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                            sheet = "porewater_data") %>%
  mutate(siteID = fct_relevel(siteID, MG.order),
         DOC = gsub("\\*", "", DOC_mg.L) %>% as.numeric(),
         SUVA = gsub("\\*", "", `SUVA_L.mg*m`) %>% as.numeric(),
         UV_254_AU = as.numeric(`UV_254_AU.cm-1`),
         sulfate_mg.L = gsub("nd", 0, sulfate_mg.L) %>% as.numeric(),
         sulfide_µg.L = gsub("nd", 0, sulfide_µg.L) %>% as.numeric()) %>%
  rename(matrixID = siteID) %>%
  select(matrixID, DOC, SUVA, UV_254_AU, sulfate_mg.L, sulfide_µg.L)


#### Data loading: RMP data ####
Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds")


#### Data management: Join RMP and porewater data ####
Hg.inc.data <- Hg.inc.data %>%
  select(coreID, matrixID, RMP_porewater) %>%
  rename(siteID = coreID)
all.data <- Hg.inc.data %>%
  left_join(porewater.data)


#### Plot: RMP of porewaters along sulfate gradient ####
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


#### Plot: RMPmatrix vs. sulfide ####
# Sulfide first
RMP.sulfide <- all.data %>%
  ggplot(aes(x = sulfide_µg.L,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Sulfide (µg/L)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = c(0.6, 0.35))
RMP.sulfide
# Sulfide spans too many orders of magnitude to accurately visualize this.
# Let's visualize this on a log scale
# To do that, let's convert the 0s to 1
all.data$sulfide_µg.L[which(all.data$sulfide_µg.L == 0)] <- 1
RMP.sulfide.log <- all.data %>%
  ggplot(aes(x = sulfide_µg.L,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_x_continuous(limits = c(-1, 4000),
                     trans = 'log10') +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Sulfide (µg/L)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.5, 0.45))
RMP.sulfide.log


#### Stats: linear model of RMP porewater vs. sulfide ####
PW.linear.model <- lm((RMP_porewater) ~ log(sulfide_µg.L, 10),
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
# Not even close to normal

shapiro.test(PW.linear.model$residuals)
# Not normally distributed, statistically speaking!

summary(PW.linear.model)



#### Plot: RMPmatrix vs sulfide with linear model ####
rSquared <- round(summary(PW.linear.model)$adj.r.squared, 2)
p.value <- with(summary(PW.linear.model), pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=F))
RMP.sulfide.plot.with.regression <- all.data %>%
  ggplot(aes(x = sulfide_µg.L,
             y = RMP_porewater,
             color = matrixID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  scale_fill_manual("black") +
  labs(x = "Sulfide (µg/L)",
       y = "Porewater RMP (%)",
       title = element_blank()) +
  theme_bw() +
  scale_x_continuous(limits = c(-1, 4000),
                     trans = 'log10') +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.7, 0.6)) +
  geom_abline(slope = coef(PW.linear.model)[[2]],
              intercept = coef(PW.linear.model)[[1]]) +
  geom_label(x = 1.5, y = 82,
             label = paste("Adjusted r2 = ", round(summary(PW.linear.model)$adj.r.squared, 3), "\n",
                           "p = ", round(p.value, 3),
                           sep = ""),
             color = "black")
RMP.sulfide.plot.with.regression





#### Plot: RMPmatrix vs. sulfate ####
# sulfate first
RMP.sulfate <- all.data %>%
  ggplot(aes(x = sulfate_mg.L,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point() +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Sulfate (µg/L)") +
  theme(axis.text.y = element_text(colour="black"),
        axis.text.x = element_text(colour="black"),
        legend.position = c(0.6, 0.35))
RMP.sulfate
# Sulfate spans too many orders of magnitude to accurately visualize this.
# Let's visualize this on a log scale
# To do that, let's convert the 0s to 1
all.data$sulfate_mg.L[which(all.data$sulfate_mg.L == 0)] <- 0.25
RMP.sulfate.log <- all.data %>%
  ggplot(aes(x = sulfate_mg.L,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_x_continuous(limits = c(-1, 30),
                     trans = 'log10') +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "Sulfate (µg/L)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.5, 0.45))
RMP.sulfate.log


#### Stats: linear model of RMP porewater vs. sulfate ####
PW.linear.model <- lm((RMP_porewater) ~ log(sulfate_mg.L, 10),
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
# Not even close to normal

shapiro.test(PW.linear.model$residuals)
# Not normally distributed, statistically speaking!

summary(PW.linear.model)




#### Plot: RMPmatrix vs sulfate with linear model ####
# rSquared <- round(summary(PW.linear.model)$adj.r.squared, 2)
# p.value <- with(summary(PW.linear.model), pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=F))
# RMP.sulfate.plot.with.regression <- all.data %>%
#   ggplot(aes(x = sulfate_mg.L,
#              y = RMP_porewater,
#              color = matrixID)) +
#   geom_smooth(method = lm ,
#               color = "black",
#               fill = "grey75",
#               se = TRUE,
#               level = 0.98) +
#   geom_point(size = 3,
#              aes(shape = matrixID)) +
#   scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
#   scale_color_manual(values = color.vector, name = "Porewater\nsource") +
#   scale_fill_manual("black") +
#   labs(x = "Sulfate (mg/L)",
#        y = "Porewater RMP (%)",
#        title = element_blank()) +
#   theme_bw() +
#   scale_x_continuous(limits = c(-1, 30),
#                      trans = 'log10') +
#   theme(axis.text.y = element_text(color = "black"),
#         axis.text.x = element_text(color = "black"),
#         legend.position = c(0.7, 0.6)) +
#   geom_abline(slope = coef(PW.linear.model)[[2]],
#               intercept = coef(PW.linear.model)[[1]]) +
#   geom_label(x = 1.5, y = 82,
#              label = paste("Adjusted r2 = ", round(summary(PW.linear.model)$adj.r.squared, 3), "\n",
#                            "p = ", round(p.value, 3),
#                            sep = ""),
#              color = "black")
# RMP.sulfate.plot.with.regression



#### Plot: RMPmatrix vs. SUVA ####
all.data %>%
  filter(!is.na(SUVA)) %>%
  ggplot(aes(x = SUVA,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "SUVA (L/mg/m)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")



#### Stats: linear model of RMP porewater vs. SUVA ####
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


#### Plot: RMPmatrix vs SUVA with linear model ####
rSquared <- round(summary(PW.linear.model)$adj.r.squared, 2)
RMP.SUVA.plot <- all.data %>%
  filter(!is.na(SUVA)) %>%
  ggplot(aes(x = SUVA,
             y = RMP_porewater,
             colour = matrixID)) +
  geom_jitter(size = 3,
             aes(shape = matrixID),
             width = 0.012) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  scale_fill_manual("black") +
  labs(x = "SUVA (L/mg/m)",
       y = "Porewater RMP (%)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.82, 0.33)) +
  geom_abline(slope = coef(PW.linear.model)[[2]],
              intercept = coef(PW.linear.model)[[1]]) +
  geom_label(x = 4, y = 82,
             label = paste("Adjusted r2 = ", round(summary(PW.linear.model)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
RMP.SUVA.plot


#### Save plot: of porewater RMP vs SUVA ####
pdf("results/incubations/RMP_porewater_SUVA.pdf",
    height = 5.5,
    width = 7)
RMP.SUVA.plot
dev.off()
# As RDS
saveRDS(object = RMP.SUVA.plot,
        file = "results/incubations/RMP_porewater_SUVA.rds")






#### Plot: RMPmatrix vs. DOC ####
all.data %>%
  filter(!is.na(DOC)) %>%
  ggplot(aes(x = DOC,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "DOC (mg/L)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")


#### Stats: linear model of RMP porewater vs. DOC ####
doc.linear.model <- lm(RMP_porewater ~ DOC,
                       data = all.data)
# Check residuals
par(mfrow = c(1,2))
plot(density(doc.linear.model$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(doc.linear.model$residuals)
qqline(doc.linear.model$residuals)
shapiro.test(doc.linear.model$residuals)
# Normally distributed, statistically speaking!
summary(doc.linear.model)
# Major effect of DOC on RMP
# Adjusted R-squared: 0.4045
# p-value: 8.458-14


#### Plot: RMPmatrix vs DOC with linear model ####
RMP.DOC.plot <- all.data %>%
  filter(!is.na(DOC)) %>%
  ggplot(aes(x = DOC,
             y = RMP_porewater,
             colour = matrixID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  scale_fill_manual("black") +
  labs(x = "DOC (mg/l)",
       y = "Porewater RMP (%)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.82, 0.33)) +
  geom_abline(slope = coef(doc.linear.model)[[2]],
              intercept = coef(doc.linear.model)[[1]]) +
  geom_label(x = 4, y = 82,
             label = paste("Adjusted r2 = ", round(summary(doc.linear.model)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
RMP.DOC.plot











#### Plot: RMPmatrix vs. UV254 absorbance ####
all.data %>%
  mutate(UV_254_AU = DOC * SUVA / 100) %>%
  filter(!is.na(UV_254_AU)) %>%
  ggplot(aes(x = UV_254_AU,
             y = RMP_porewater,
             group = matrixID,
             color = matrixID)) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(name = "Spiking matrix",
                     values = color.vector) +
  theme_bw() +
  labs(y = "Porewater RMP",
       x = "UV 254 (cm-1)") +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = "none")

#### Stats: linear model of RMPmatrix vs. UV254 absorbance ####
UV.linear.model <- lm(RMP_porewater ~ UV_254_AU,
                      data = all.data)
# Check residuals
par(mfrow = c(1,2))
plot(density(UV.linear.model$residuals),
     main="Density plot of residuals",
     ylab="Density",
     xlab="Residuals")
# QQ-normal plot
qqnorm(UV.linear.model$residuals)
qqline(UV.linear.model$residuals)
shapiro.test(UV.linear.model$residuals)
# Normally distributed, statistically speaking!
summary(UV.linear.model)
# Major effect of DOC on RMP
# Adjusted R-squared: 0.3755
# p-value: 3.504-10


#### Plot: RMPmatrix vs UV254 with linear model ####
RMP.UV.plot <- all.data %>%
  mutate(UV_254_AU = DOC * SUVA / 100) %>%
  filter(!is.na(UV_254_AU)) %>%
  ggplot(aes(x = UV_254_AU,
             y = RMP_porewater,
             colour = matrixID)) +
  geom_smooth(method = lm ,
              color = "black",
              fill = "grey75",
              se = TRUE,
              level = 0.98) +
  geom_point(size = 3,
             aes(shape = matrixID)) +
  scale_shape_manual(values = point.vector, name = "Porewater\nsource") +
  scale_color_manual(values = color.vector, name = "Porewater\nsource") +
  scale_fill_manual("black") +
  labs(x = "UV254 (cm-1)",
       y = "Porewater RMP (%)",
       title = element_blank()) +
  theme_bw() +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = c(0.3, 0.75)) +
  geom_abline(slope = coef(UV.linear.model)[[2]],
              intercept = coef(UV.linear.model)[[1]]) +
  geom_label(x = 1.2, y = 75,
             label = paste("Adjusted R2 = ", round(summary(UV.linear.model)$adj.r.squared, 3), "\n",
                           "p << 0.0001",
                           sep = ""),
             color = "black")
RMP.UV.plot







#### Plot em all out ####
(RMP.sulfide.log + RMP.SUVA) / (RMP.DOC.plot + RMP.UV)



#### Save out the ones with no correlation ####
saveRDS(object = RMP.sulfide.plot.with.regression,
        file = "results/incubations/RMP_porewater_sulfide.rds")
saveRDS(object = RMP.sulfate.log,
        file = "results/incubations/RMP_porewater_sulfate.rds")
saveRDS(object = RMP.DOC.plot,
        file = "results/incubations/RMP_porewater_doc.rds")
saveRDS(object = RMP.UV.plot,
        file = "results/incubations/RMP_porewater_uv.rds")

