#### code/readBasedAnalysis/graftm_phyla_plots.R ####
# Benjamin D. Peterson


#### Clean up first ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### Read in metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Replicate metadata
rep.metadata.vector <- metadata$replicate
names(rep.metadata.vector) <- metadata$metagenomeID



#### Read in data ####
tax.data <- readRDS("dataEdited/readBased_analysis/16S/graftM/taxonomy_counts.rds") %>%
  filter(grepl("KMBP005", metagenomeID))



#### Clean up data ####
phyla.data <- tax.data %>%
  group_by(metagenomeID, phylum) %>%
  summarize(rel.abundance = sum(rel.abundance,
                                na.rm = TRUE)) %>%
  filter(!is.na(phylum)) %>%
  mutate(siteID = site.metadata.vector[metagenomeID]) %>%
  ungroup()


#### Identify abundant phyla ####
needed.phyla.to.plot <- phyla.data %>%
  filter(rel.abundance >= 0.05) %>%
  select(phylum) %>%
  unlist(use.names = FALSE) %>%
  unique()


#### Make color vector ####
CB.color.vector <- c("gray75",
                     readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds"),
                     "gray", "gray35")
phylum.colors <- data.frame(phylum = needed.phyla.to.plot %>% sort(),
                            phylum.color = CB.color.vector)

#### Generate plot ####
pdf("results/figures/S10/graftm_phyla.pdf",
    width = 6,
    height = 3)
phyla.data %>%
  filter(phylum %in% needed.phyla.to.plot) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = rep.metadata.vector[metagenomeID],
             y = rel.abundance,
             fill = phylum)) +
  geom_col() +
  scale_fill_manual(values = phylum.colors$phylum.color,
                    name="Phylum",
                    breaks = phylum.colors$phylum) +
  facet_grid(~siteID) +
  theme_bw() +
  xlab("") +
  ylab("Relative abundance") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 1))
dev.off()
 