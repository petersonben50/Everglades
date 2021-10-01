#### Clean up first ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)



#### Read in metadata ####

metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")

# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Replicate metadata
rep.metadata.vector <- metadata$replicate
names(rep.metadata.vector) <- metadata$metagenomeID





#### Generate vector of correct order of samples ####

MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")



#### Read in data ####

tax.data <- readRDS("dataEdited/2019_readBasedAnalysis/16S/graftM/taxonomy_counts.rds") %>%
  filter(grepl("KMBP005", metagenomeID))


phyla.data <- tax.data %>%
  group_by(metagenomeID, phylum) %>%
  summarize(rel.abundance = sum(rel.abundance,
                                na.rm = TRUE)) %>%
  filter(!is.na(phylum)) %>%
  mutate(siteID = site.metadata.vector[metagenomeID]) %>%
  ungroup()
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
pdf("results/2019_readBasedAnalysis/graftm_phyla.pdf",
    width = 8,
    height = 4)
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
 