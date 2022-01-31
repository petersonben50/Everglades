#### code/assembly_analysis/hgcA/hgcA_graftM_ratios.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(patchwork)
library(readxl)
library(tidyverse)
# source("code/setup_PW_core_order_color_points.R")


##### Read in MG metadata ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID

medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID


#### Read in hgcA data ####
metagenome.ids.to.include <- names(medium.renaming.vector)[grepl("KMBP005", names(medium.renaming.vector))]
hgcA.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                      stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_class")) %>%
  arrange(classification) %>%
  filter(seqID %in% readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt")) %>%
  select(clusterName, all_of(metagenome.ids.to.include)) %>%
  gather(key = metagenomeID,
         value = abundance,
         -1) %>%
  group_by(metagenomeID, clusterName) %>%
  summarise(hgcA.abundance = sum(abundance))


#### Read in taxonomic composition data ####
tax.data <- readRDS("dataEdited/readBased_analysis/16S/graftM/taxonomy_counts.rds") %>%
  filter(grepl("KMBP005", metagenomeID))

# Pull out phyla of interest
tax.data.phyla <- tax.data %>%
  filter(phylum %in% c("Chloroflexi")) %>%
  rename(clusterName = phylum) %>%
  group_by(metagenomeID, clusterName) %>%
  summarise(tax.abundance = sum(rel.abundance,
                                na.rm = TRUE)*100)

# Pull out classes of interest
tax.data.classes <- tax.data %>%
  filter(class %in% c("Aminicenantia")) %>%
  mutate(class = "Aminicenantes") %>%
  rename(clusterName = class) %>%
  group_by(metagenomeID, clusterName) %>%
  summarise(tax.abundance = sum(rel.abundance,
                                na.rm = TRUE)*100)

# Pull out orders of interest
tax.data.orders <- tax.data %>%
  filter(order %in% c("Syntrophobacterales")) %>%
  rename(clusterName = order) %>%
  group_by(metagenomeID, clusterName) %>%
  summarise(tax.abundance = sum(rel.abundance,
                                na.rm = TRUE)*100)

# Combine tax data
tax.data.of.interest <- rbind(tax.data.phyla, tax.data.classes, tax.data.orders)
rm(tax.data.phyla, tax.data.classes, tax.data.orders)


#### Combine data ####

all.data.of.interest <- hgcA.data %>%
  right_join(tax.data.of.interest) %>%
  mutate(percent_methylators = hgcA.abundance  / tax.abundance * 100)
  