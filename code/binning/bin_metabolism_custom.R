#### code/2019_binning/bin_metabolism_custom.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(readxl)
library(tidyverse)


#### Read in taxa data ####
tax.data <- read.table("dataEdited/2019_binning/binning_initial/binsGood/taxonomy_summary.txt",
                       sep = '\t',
                       col.names = c("binID", "taxonomy"))


#### Read in MHC data ####
MHC.data <- read.table("dataEdited/2019_binning/MHCs/heme_count_bins.tsv",
                       sep = '\t',
                       header = TRUE)
MHC.data.grouped <- MHC.data %>%
  group_by(binID) %>%
  summarize(mhcGeneCount = n())
all.data <- full_join(tax.data,
                      MHC.data.grouped)

#### Read in metabolism data ####
metabolism.data <- read.table("dataEdited/2019_binning/batch_HMM_output/bin_counts/all_bin_hits.tsv",
                              sep = '\t',
                              header = TRUE) %>%
  group_by(binID, proteinName) %>% 
  mutate(geneID = paste0(geneID, collapse = ",")) %>%
  slice(1) %>%
  spread(key = proteinName,
         value = geneID)
all.data <- all.data %>%
  right_join(metabolism.data)


#### Save out this data frame ####
write.csv(metabolism.data,
          "dataEdited/2019_binning/batch_HMM_output.csv",
          row.names = FALSE)
