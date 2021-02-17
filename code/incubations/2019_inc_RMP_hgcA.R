#### code/incubations/2019_inc_RMP_hgcA.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")


#### Read in data ####
RMP.data <- readRDS("dataEdited/2019_incubations/rel_methylation_microbes.rds")
hgcA.data <- readRDS("dataEdited/2019_analysis_assembly/hgcA/hgcA_abundance_site.rds")


#### Prep color and point vector ####
color.vector <- cb.translator[1:6]
names(color.vector) <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")
point.vector <- c(9, 18:15, 25)
names(point.vector) <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")


#### Prep data for plotting ####
plotting.data.hgcA <- hgcA.data %>%
  group_by(siteID) %>%
  summarise(coverage_mean = mean(coverage),
            coverage_sd = sd(coverage),
            coverage_count = n(),
            coverage_se = coverage_sd / sqrt(coverage_count))
plotting.data.inc <- RMP.data %>%
  rename(siteID = coreID) %>%
  group_by(siteID) %>%
  summarise(rel_meth_spike_mean = mean(rel_meth_spike * 100),
            rel_meth_spike_sd = sd(rel_meth_spike * 100),
            rel_meth_spike_count = n()) %>%
  mutate(rel_meth_spike_se = rel_meth_spike_sd / sqrt(rel_meth_spike_count)) %>%
  ungroup()
plotting.data <- full_join(plotting.data.hgcA,
                           plotting.data.inc)
rm(plotting.data.inc,
   plotting.data.hgcA)


#### Plot abundance of hgcA in sediment against relative methylation ####
pdf("results/2019_analysis_assembly/relative_methylation_vs_hgcA_coverage_sediment.pdf",
    width = 6,
    height = 4)
plotting.data %>%
  ggplot(aes(x = coverage_mean,
             y = rel_meth_spike_mean,
             colour = siteID)) +
  geom_errorbar(aes(ymin = rel_meth_spike_mean - rel_meth_spike_se,
                    ymax = rel_meth_spike_mean + rel_meth_spike_se),
                colour = "black",
                width = 0.25) +
  geom_errorbarh(aes(xmin = coverage_mean - coverage_se,
                     xmax = coverage_mean + coverage_se),
                 colour = "black",
                 height = 2.5) +
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
        axis.text.x = element_text(color = "black"))
dev.off()