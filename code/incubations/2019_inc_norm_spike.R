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

PW.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8", "FI", "CYSTEINE", "CONTROL")
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Prep color vector ####
color.vector.core <- cb.translator[1:6]
names(color.vector.core) <- MG.order
color.vector.PW <- c(cb.translator[1:6], cb.translator["reddishpurple"], "grey50", "grey75")
names(color.vector.PW) <- PW.order


#### Read in incubation data and order it ####
inc.Hg.data <- read_xlsx("dataRaw/geochem/2019/December 2019 field trip_synthesis of core experiments_v2.xlsx",
                         sheet = "incubation_results_cleaned") %>%
  mutate(coreID = fct_relevel(coreID, MG.order)) %>%
  mutate(matrixID = fct_relevel(matrixID, PW.order))



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
  scale_fill_manual(values = color.vector.PW) +
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
  mutate(rel_meth_spike = (MeHg_fract_spike) / (max_meth)) %>%
  mutate(matrixID = fct_relevel(matrixID, PW.order))


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
  scale_fill_manual(values = color.vector.PW) +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(y = "Methylated spike (%)") +
  theme(axis.text.y = element_text(colour="black"))
# dev.off()


#### Plot RMP as points with lines connecting them ####
# pdf("results/2019_incubations/spikeNormalized_all.pdf",
#     width = 6,
#     height = 4)
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
  scale_color_manual(values = color.vector.PW) +
  theme_classic() +
  labs(y = "RMP",
       x = "Sediment core site",
       title = "Relative methylation potential (RMP)\nof porewater solutions") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
# dev.off()



#### Read in geochemical data for comparison ####
PW.geochem.data <- read_xlsx("dataEdited/geochem/geochem_data_2019.xlsx",
                             sheet = "PW_SW_geochem") %>%
  filter(medium == "PW") %>%
  select(siteID, SUVA_254_L.mg.m, Sulfide_µg.L) %>%
  rename(matrixID = siteID) %>%
  filter(SUVA_254_L.mg.m != "NA")
all.data.suva <- PW.geochem.data %>%
  left_join(all.data %>% mutate(matrixID = as.character(matrixID))) %>%
  mutate(suva = as.numeric(SUVA_254_L.mg.m))

# Add in WCA data
all.data.suva$siteLocation <- "WCA"
all.data.suva$siteLocation[grep("LOX", all.data.suva$matrixID)] <- "LOX"

#### Plot RMP against the SUVA content ####
pdf("results/2019_incubations/spikeNormalized_suva.pdf",
    width = 6,
    height = 4)
all.data.suva %>%
  group_by(coreID, matrixID, suva, siteLocation) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100) %>%
  ggplot(aes(x = suva,
             y = meth_spike_per_mean,
             group = coreID)) + 
  geom_point(aes(color = coreID,
                 shape = siteLocation)) +
  geom_line(aes(color = coreID)) +
  scale_color_manual(values = color.vector.PW) +
  theme_classic() +
  labs(y = "RMP",
       x = "SUVA (L/mg/m)",
       title = "Relative methylation potential (RMP)\nof sediment cores relative to SUVA content") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
dev.off()

#### Plot RMP against the sulfide content ####
all.data.suva %>%
  group_by(coreID, matrixID, Sulfide_µg.L, siteLocation) %>%
  summarise(meth_spike_per_mean = mean(rel_meth_spike) * 100) %>%
  ggplot(aes(x = log(Sulfide_µg.L, 10),
             y = meth_spike_per_mean,
             group = coreID)) + 
  geom_point(aes(color = coreID,
                 shape = siteLocation)) +
  geom_line(aes(color = coreID)) +
  scale_color_manual(values = color.vector.PW) +
  theme_classic() +
  labs(y = "RMP",
       x = "SUVA (L/mg/m)",
       title = "Relative methylation potential (RMP)\nof sediment cores relative to SUVA content") +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
