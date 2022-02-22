#### code/assembly_analysis/merB/abundance_merB.R ####
# Benjamin D. Peterson

# Scripts to look at abundance of merB in a myriad of ways
# Can't believe I wrote "a myriad of ways"

#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")


#### List of merB sequences ####
merB.list <- readLines("dataEdited/assembly_analysis/merB/merB_raw.txt")
names(merB.list) <- paste(strsplit(merB.list, "_") %>% sapply("[", 1),
                          strsplit(merB.list, "_") %>% sapply("[", 2),
                          sep = "_")

#### Read in abundance and phylogenetic info ####
all.data <- read.csv("dataEdited/assembly_analysis/merB/depth/merB_coverage_raw_scgNormalization.csv",
                     stringsAsFactors = FALSE) %>%
  mutate(seqID = merB.list[scaffoldID]) %>%
  filter(seqID %in% merB.list)


#### Read in metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID

metagenomes.names.to.use <- grep("KMBP00",
                                 metadata$metagenomeID,
                                 value = TRUE)



#### Prep color vector ####
color.vector <- color.vector[1:6]


#### Plot abundance of merB ####
merB.plot <- all.data %>%
  gather(key = MG,
         value = coverage,
         c(metagenomes.names.to.use)) %>%
  group_by(MG) %>%
  summarise(coverage = sum(coverage)) %>%
  mutate(siteID = site.metadata.vector[MG],
         medium = medium.metadata.vector[MG]) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  filter(medium == "sediment") %>%
  ggplot(aes(x = siteID,
             y = coverage,
             color = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  ylab("merB coverage") +
  xlab("Peat core ID") +
  theme_bw()


#### Save out merB plot ####
pdf("results/figures/S13_merB_abundance.pdf",
    height = 3,
    width = 5)
merB.plot
dev.off()
