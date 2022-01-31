#### code/assembly_analysis/methanogenesis/methanogen_abundance.R ####
# Benjamin D. Peterson



#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(ggplot2)
library(patchwork)
library(readxl)
library(tidyverse)
source("code/setup_PW_core_order_color_points.R")



##### Read in MG metadata ####

MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID

medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID



##### Read in depth info ####

depth.data <- read.csv("dataEdited/assembly_analysis/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.renaming.vector[sampleID]) %>%
  mutate(medium = medium.renaming.vector[sampleID])




#### Read in mcrA list ####

mcrA.scaffold.list <- readLines("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)


#### Plot abundance with taxonomy ####
# Color vector first
color.vector.mcrA <- c(cb.translator["reddishpurple"],
                       cb.translator["orange"],
                       cb.translator["vermillion"])
names(color.vector.mcrA) <- read_xlsx("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_clusters.xlsx") %>%
  select(clusterName) %>%
  unlist(use.names = FALSE) %>%
  unique()

# Prep data for plotting
depth.tax <- depth.data %>%
  filter(scaffoldID %in% mcrA.scaffold.list) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_clusters.xlsx") %>%
              mutate(scaffoldID = seqID %>% strsplit("_[1-9]+") %>% sapply("[", 1))) %>%
  group_by(sampleID, siteID, medium, clusterName) %>%
  summarize(coverage = round(sum(coverage), 4)) %>%
  ungroup() %>%
  filter(medium == "sediment") %>%
  group_by(siteID, clusterName) %>%
  summarize(coverage = mean(coverage)) %>%
  ungroup() %>%
  mutate(siteID = fct_relevel(siteID, MG.order))

# Generate plot
mcrA.plot <- depth.tax %>%
  ggplot(aes(x = siteID,
             y = coverage,
             fill = clusterName)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color.vector.mcrA,
                    name = "Taxonomic group") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Gene coverage\n(per 100X coverage of SCGs)") +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        legend.position = c(0.3, 0.7))

# Save out plot
saveRDS(mcrA.plot,
        "results/metagenomes/assembly/methanogenesis/mcrA_abundance.rds")






#### Read in marker lists ####

list.of.marker.lists <- list.files(path = "dataEdited/assembly_analysis/metabolicProteins/methanogenesis/markers/",
                                   pattern = "methan_mark")
for (marker.of.interest in 1:length(list.of.marker.lists)) {

  scaffolds.of.interest = readLines(paste("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/markers/",
                                          list.of.marker.lists[marker.of.interest],
                                          sep = "")) %>%
    strsplit("_[1-9]+") %>%
    sapply("[", 1)
  marker.df.temp <- data.frame(scaffoldID = scaffolds.of.interest,
                               geneName = rep(paste("marker_", marker.of.interest,
                                                    sep = ""),
                                              length(scaffolds.of.interest)))

  if (marker.of.interest == 1) {
    marker.df <- marker.df.temp
  } else {
    marker.df <- rbind(marker.df,
                       marker.df.temp)
  }
  rm(marker.df.temp)
}

marker.depth <- left_join(marker.df,
                          depth.data)


#### Plot data ####
marker.depth %>%
  filter(medium == "sediment") %>%
  group_by(scaffoldID, geneName, siteID) %>%
  summarise(coverage = mean(coverage)) %>%
  ungroup() %>%
  group_by(geneName, siteID) %>%
  summarise(coverage = sum(coverage)) %>%
  ungroup() %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = coverage)) +
  geom_line(mapping = aes(group = geneName)) +
  geom_point() +
  theme_classic() +
  stat_summary(geom = "point", fun = "mean",
               col = "black", fill = "red",
               size = 3, shape = 24) +
  labs(title = "Coverage of methanogenic marker genes",
       y = "Normalized to SCG coverage (100X)",
       x = "Site ID")
