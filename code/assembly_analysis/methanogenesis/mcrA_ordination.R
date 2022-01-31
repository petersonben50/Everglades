#### code/assembly_analysis/sulfur/dsrA_ordination.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggpubr)
library(phyloseq)
library(readxl)
library(tidyverse)
library(vegan)
source("code/setup_PW_core_order_color_points.R")


#### Read in metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx") %>%
  filter(grepl("KMBP005",
               x = metagenomeID))
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID


#### Read in dsrA lists ####
mcrA.scaffolds <- readLines("dataEdited/assembly_analysis/metabolicProteins/methanogenesis/mcrA/mcrA_derep_list.txt") %>%
  strsplit("_[1-9]+") %>%
  sapply("[", 1)


##### Read in depth info ####
depth.data <- read.csv("dataEdited/assembly_analysis/metabolicProteins/depth/metabolicProtein_depth_clean.csv",
                       stringsAsFactors = FALSE) %>%
  gather(key = sampleID,
         value = coverage,
         -1) %>%
  mutate(siteID = site.metadata.vector[sampleID]) %>%
  mutate(medium = medium.metadata.vector[sampleID]) %>%
  rename(MG = sampleID) %>%
  filter(scaffoldID %in% mcrA.scaffolds,
         medium == "sediment") %>%
  select(-medium) %>%
  mutate(siteID = fct_relevel(siteID, MG.order))



#### Calculate "richness" ####
richness.data <- depth.data %>%
  filter(coverage > 0) %>%
  group_by(MG, siteID) %>%
  summarise(unique_sequences = n())
richness.plot <- richness.data %>%
  ggplot(aes(x = siteID,
             y = unique_sequences,
             color = siteID)) +
  geom_point() +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  ylim(c(0, 30)) +
  ylab("Richness") +
  theme_bw() +
  theme(legend.position = "none",
        # legend.position = c(0.8, 0.4),
        axis.title.x = element_blank())
richness.plot


#### Evenness - Shannon diversity ####
total.data <- depth.data %>%
  group_by(MG) %>%
  summarise(total.coverage = sum(coverage))
evenness.data <- depth.data %>%
  left_join(total.data) %>%
  mutate(relative_coverage = coverage / total.coverage) %>%
  select(MG, scaffoldID, relative_coverage) %>%
  spread(key = scaffoldID,
         value = relative_coverage)
row.names(evenness.data) <- evenness.data$MG
evenness.data <- evenness.data %>%
  select(-MG) %>%
  as.matrix()

# Generate phyloseq object
mcrA.phyloseq <- otu_table(evenness.data,
                           taxa_are_rows = FALSE)

# Calculate Shannon diversity
diversity.vector <- diversity(mcrA.phyloseq,
                              index = "shannon")
diversity.df <- data.frame("MG" = names(diversity.vector),
                           "diversity" = diversity.vector) %>%
  mutate(siteID = site.metadata.vector[MG])
row.names(diversity.df) <- NULL

# Plot Shannon diversity
diversity.df %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = diversity)) +
  geom_point() +
  ylim(c(0, 4)) +
  theme_bw()


#### Shannon's equitability ####
E.data <- diversity.df %>%
  left_join(richness.data %>% select(-siteID)) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(E = diversity / log(unique_sequences))
E.plot <- E.data %>%
  ggplot(aes(x = siteID,
             y = E,
             color = siteID)) +
  geom_point() +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
E.plot


#### Plot em ####
ggarrange(richness.plot,
          E.plot)



#### Statistical analysis of evenness ####
richness.model <- aov(unique_sequences ~ siteID, data = richness.data)
summary(richness.model)
TukeyHSD(richness.model, conf.level=.95) 

evenness.model <- aov(E ~ siteID, data = E.data)
summary(evenness.model)
TukeyHSD(evenness.model, conf.level=.95) 





#### Prepare data for ordination ####

#make metadata a phyloseq class object
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[, "metagenomeID"] %>%
  unlist()
meta.phylo <- sample_data(metadata)

##make biom for phyloseq
all.phyloseq <- merge_phyloseq(mcrA.phyloseq, meta.phylo)


#### Generate ordination ####
mcrA.bray <- ordinate(all.phyloseq,
                      method = "PCoA",
                      distance = "bray")

bc.ord.mcrA <- plot_ordination(all.phyloseq,
                               mcrA.bray,
                               color = "siteID",
                               shape = "siteID") +
  scale_color_manual(values = color.vector[1:6],
                     name = "Site ID") +
  scale_shape_manual(values = point.vector[1:6],
                     name = "Site ID") +
  geom_point(size=5) +
  ylim(c(-0.5, 0.25)) +
  xlim(c(-0.5, 0.6)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.2, 0.3))
bc.ord.mcrA


#### Save out ordination ####
saveRDS(bc.ord.mcrA,
        "results/metagenomes/assembly/methanogenesis/mcrA_ordination.rds")
