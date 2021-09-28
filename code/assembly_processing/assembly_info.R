#### code/2019_assembly_info.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(tidyverse)



#### Read in data ####
assembly.info <- read.table("dataEdited/assemblies/reports/all_assemblies_stats.txt",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  select(assemblyID, n, L50, N50, sum)

ORF.counts <- read.table("dataEdited/assemblies/reports/ORF_counts.tsv",
                         stringsAsFactors = FALSE,
                         header = TRUE)



#### Combine data ####
all.info <- full_join(assembly.info,
                      ORF.counts) %>%
  filter(grepl("Sed99",
               assemblyID))



#### Write out data table ####
write.csv(x = all.info,
          "results/2019_assembly_info.csv",
          row.names = FALSE)
