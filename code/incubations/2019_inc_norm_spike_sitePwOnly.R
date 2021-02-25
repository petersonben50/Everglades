#### code/incubations/2019_inc_norm_spike.R ####
# Benjamin D. Peterson

#### Get set up #####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(AICcmodavg)
library(lme4)
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Generate vector of correct order of samples along sulfate gradient ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Prep color vector ####
color.vector.core <- cb.translator[1:6]
names(color.vector.core) <- MG.order


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                         sheet = "incubation_results_cleaned") %>%
  filter(matrixID %in% MG.order) %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) %>%
  mutate(matrixID = fct_relevel(matrixID, MG.order))


#### Plot methylation data by site, grouped within spiking solution ####

# pdf("results/2019_incubations/spike_faceted.pdf",
#     width = 9,
#     height = 6)
inc.Hg.data %>%
  mutate(coreID = paste("Core source: ",
                        coreID,
                        sep = "")) %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(MeHg_fract_spike) * 100,
            meth_spike_per_sd = sd(MeHg_fract_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(coreID = fct_relevel(coreID, paste("Core source: ",
                                            MG.order,
                                            sep = ""))) %>%
  ggplot(aes(x = matrixID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = matrixID)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~coreID, nrow = 2) +
  scale_fill_manual(values = color.vector.core) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))
# dev.off()



#### Plot methylation data, normalizing to highest % within the spiking solution ####

inc.Hg.data.sum <- inc.Hg.data %>%
  group_by(coreID) %>%
  summarise(max_meth = max(MeHg_fract_spike))
all.data <- full_join(inc.Hg.data,
                      inc.Hg.data.sum) %>%
  mutate(rel_meth_spike = (MeHg_fract_spike) / (max_meth))



# pdf("results/2019_incubations/spikeNormalized_faceted.pdf",
#     width = 9,
#     height = 6)
all.data %>%
  mutate(coreID = paste("Core source: ",
                        coreID,
                        sep = "")) %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100,
            meth_spike_per_sd = sd(rel_meth_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ungroup() %>%
  mutate(coreID = fct_relevel(coreID, paste("Core source: ",
                                            MG.order,
                                            sep = ""))) %>%
  ggplot(aes(x = matrixID,
             y = meth_spike_per_mean,
             width = 0.8,
             fill = matrixID)) + 
  geom_bar(position="stack", stat="identity") +
  geom_errorbar(aes(ymin = meth_spike_per_mean - meth_spike_per_se,
                    ymax = meth_spike_per_mean + meth_spike_per_se),
                colour = "black",
                width = 0.33) +
  facet_wrap(~coreID, nrow = 2) +
  scale_fill_manual(values = color.vector.core) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))
# dev.off()


#### Plot RMP as points with lines connecting them ####
pdf("results/2019_incubations/spikeNormalized_RMP_sitePwOnly.pdf",
    width = 6,
    height = 4)
all.data %>%
  group_by(coreID, matrixID) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100,
            meth_spike_per_sd = sd(rel_meth_spike),
            meth_spike_per_count = n(),
            meth_spike_per_se = meth_spike_per_sd / sqrt(meth_spike_per_count) * 100) %>%
  ggplot(aes(x = matrixID,
             y = meth_spike_per_mean,
             group = coreID)) + 
  geom_point(aes(color = coreID)) +
  geom_line(aes(color = coreID)) +
  scale_color_manual(values = color.vector.core) +
  theme_classic() +
  labs(y = "RMP",
       x = "Porewater origin",
       title = "Relative methylation potential (RMP) of porewaters") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()


#### Save out normalized data ####
all.data %>%
  saveRDS("dataEdited/2019_incubations/rel_methylation_PW.rds")
