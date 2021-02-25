#### code/incubations/2019_inc_norm_core.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(AICcmodavg)
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
PW.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8", "FI", "CYSTEINE", "CONTROL")


#### Prep color vector ####
color.vector <- cb.translator[1:6]
names(color.vector) <- MG.order
color.vector.PW <- c(cb.translator[1:6], cb.translator["reddishpurple"], "grey50", "grey75")
names(color.vector.PW) <- PW.order



#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                         sheet = "incubation_results_cleaned") %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) %>%
  mutate(matrixID = fct_relevel(matrixID, PW.order))


#### Plot methylation data by site, grouped within spiking solution ####

pdf("results/2019_incubations/coreFaceted_MeHgProduction.pdf",
    width = 9,
    height = 6)
inc.Hg.data %>%
  mutate(matrixID = paste("Spiking matrix\nsource: ",
                           matrixID,
                           sep = "")) %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(MeHg_fract_spike) * 100,
            meth_spike_per_sd = sd(MeHg_fract_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) %>%
  mutate(matrixID = fct_relevel(matrixID, paste("Spiking matrix\nsource: ",
                                                PW.order,
                                                sep = ""))) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~matrixID, nrow = 2) +
  labs(y = "Methylated spike (%)") +
  scale_fill_manual(values = color.vector) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(color = "black"))
dev.off()






#### Plot methylation data, normalizing to highest % within the spiking solution ####

inc.Hg.data.sum <- inc.Hg.data %>%
  group_by(matrixID) %>%
  summarise(max_meth = max(MeHg_fract_spike))
all.data <- full_join(inc.Hg.data,
                      inc.Hg.data.sum) %>%
  mutate(rel_meth_spike = (MeHg_fract_spike) / (max_meth))


pdf("results/2019_incubations/coreFaceted_MeHgProductionNormalized.pdf",
    width = 9,
    height = 6)
all.data %>%
  mutate(matrixID = paste("Spiking matrix\nsource: ",
                          matrixID,
                          sep = "")) %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100,
            meth_spike_per_sd = sd(rel_meth_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) %>%
  mutate(matrixID = fct_relevel(matrixID, paste("Spiking matrix\nsource: ",
                                                PW.order,
                                                sep = ""))) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = coreID)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~matrixID, nrow = 2) +
  scale_fill_manual(values = color.vector) +
  labs(y = "Methylated spike (%)") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()



#### Plot RMP as individual points connected within PW source ####
pdf("results/2019_incubations/coreNormalized_RMP_all.pdf",
    width = 6,
    height = 4)
all.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100,
            meth_spike_per_sd = sd(rel_meth_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             group = matrixID)) + 
  geom_point(aes(color = matrixID)) +
  geom_line(aes(color = matrixID)) +
  scale_color_manual(values = color.vector.PW) +
  labs(y = "RMP",
       x = "Sediment core site",
       title = "Relative methylation potential (RMP)\nof sediment cores") +
  theme_classic() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()

# Do with only site water
pdf("results/2019_incubations/coreNormalized_RMP_siteWater.pdf",
    width = 6,
    height = 4)
all.data %>%
  filter(matrixID %in% MG.order) %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100,
            meth_spike_per_sd = sd(rel_meth_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ggplot(aes(x = coreID,
             y = meth_spike_per_mean,
             group = matrixID)) + 
  geom_point(aes(color = matrixID)) +
  geom_line(aes(color = matrixID)) +
  scale_color_manual(values = color.vector) +
  theme_classic() +
  labs(y = "RMP",
       x = "Sediment core site",
       title = "Relative methylation potential (RMP)\nof sediment cores with in situ waters") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()


#### Save out this dataframe ####
all.data %>%
  saveRDS("dataEdited/2019_incubations/rel_methylation_microbes.rds")
