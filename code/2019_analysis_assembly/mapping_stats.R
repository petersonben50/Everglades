#### code/2019_analysis_assembly/mapping_stats.R ####
# Benjamin D. Peterson
# Analysis in 2019_analysis_assembly/mapping_assembly.md

#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(readxl)
library(tidyverse)


#### Read in data ####

# Read in overall mapping of each metagenome
uniq.mapped.reads <- read.table("dataEdited/2019_analysis_assembly/reports/uniq_mapped_reads_MG.tsv",
                                sep = "\t",
                                col.names = c("metagenomeID", "forwardReads",
                                              "reverseReads", "singleReads"))

# Read in metagenome read counts
all.reads <- read.table("dataEdited/metagenomes/reports/metagenome_read_count.tsv",
                        sep = "\t",
                        header = TRUE) %>%
  mutate(totalReads = forwardReads + reverseReads + singleReads + mergedReads) %>%
  select(metagenomeID, totalReads)

# Read in metadata
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")

#### Combine data and calculate fraction mapping to all  ####

all.data <- MG.metadata %>%
  right_join(uniq.mapped.reads) %>%
  left_join(all.reads) %>%
  mutate(percentMappingAllMGs = (forwardReads + reverseReads + singleReads) / totalReads * 100) %>%
  select(-sampleID)
