#### code/assembly_analysis/hgcA/aggregate_hgcA_information.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(readxl)
library(tidyverse)



#### Make hgcA table ####
hgcA.list <- readLines("dataEdited/assembly_analysis/hgcA/hgcA.txt")
hgcA.true.list <- readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt")
hgcA.df <- data.frame(seqID = hgcA.list,
                      scaffoldID = paste(hgcA.list %>% strsplit("_") %>% sapply("[", 1),
                                         hgcA.list %>% strsplit("_") %>% sapply("[", 2),
                                         sep = "_"))

#### Prep pplacer classification ####

taxonomy.table <- read.csv("dataEdited/assembly_analysis/hgcA/classification/hgcA_taxonomy_table.csv",
                           stringsAsFactors = FALSE) %>%
  mutate(pplacer_classification = paste(phylum, class, order,
                                        family, genus, species,
                                        sep = ";") %>% 
           strsplit(";Unclassified") %>% sapply("[", 1)) %>%
  select(seqID, pplacer_classification)



#### Manual classification ####

manual.taxonomy <- read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                             sheet = "seq_class") %>%
  rename(manual_classification = classification) %>%
  select(-c(notes, clusterName))



#### Prep site information

metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID



#### Read in coverage data ####
depth.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                       stringsAsFactors = FALSE) %>%
  select(-c(scaffoldID, length, total)) %>%
  gather(key = metagenomeID,
         value = coverage,
         -1) %>%
  filter(grepl("KMBP005", metagenomeID)) %>%
  mutate(metagenomeID = paste(site.metadata.vector[metagenomeID],
                              ";", metagenomeID,
                              sep = "")) %>%
  spread(key = metagenomeID,
         value = coverage)



#### Prep binning information ####
binning.data <- read.table("/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/bin_analysis/binning_initial/binsFinal/hgcA_bin_key.tsv",
                           sep = '\t',
                           header = TRUE) %>%
  select(binID, hgcA) %>%
  left_join(readRDS("dataEdited/assembly_analysis/hgcA/dereplication/hgcA_derep_key.rds")) %>%
  left_join(read_xlsx("dataEdited/bin_analysis/binning_initial/binsGood/taxonomy_summary.xlsx")) %>%
  rename(seqID = hgcA_rep) %>%
  select(-hgcA) %>%
  rename(binTaxonomy = taxID)


#### Prep hgcB information ####
# hgcB.data <- readLines("dataEdited/assembly_analysis/hgcA/hgcB/hgcB_good.txt")
hgcB.data <- read_xlsx("dataEdited/assembly_analysis/hgcA/hgcB/hgcB_notes.xlsx")
hgcB.df <- hgcB.data %>%
  mutate(scaffoldID = paste(hgcB_ID %>% strsplit("_") %>% sapply("[", 1),
                            hgcB_ID %>% strsplit("_") %>% sapply("[", 2),
                            sep = "_")) %>%
  inner_join(hgcA.df) %>%
  select(-scaffoldID)
rm(hgcB.data)



#### Prep paralog information ####
paralog.df <- data.frame(seqID = hgcA.list,
                         isParalog = !(hgcA.list %in% hgcA.true.list))


#### Aggregate all data ####
hgcA.data <- paralog.df %>%
  left_join(manual.taxonomy) %>%
  left_join(taxonomy.table) %>%
  left_join(binning.data) %>%
  left_join(hgcB.df) %>%
  left_join(depth.data)


#### Save data out ####
write.csv(hgcA.data,
          "dataEdited/assembly_analysis/hgcA/hgcA_dataTable.csv",
          row.names = FALSE)
