#### code/2019_analysis_assembly/clean_hgcA_coverage.R ####
# Written by Benjamin D. Peterson



#### Get ready. Get set ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(stringr)
library(tidyverse)



#### Read in normalization vector ####

normalized.coverage.vector <- readRDS("dataEdited/metagenomes/reports/metagenome_normalization_vector.rds")



#### List of file names with depth info ####

list.o.depth.files <- list.files(path = "dataEdited/2018_analysis_assembly/hgcA/depth",
                                 pattern = "hgcA_depth",
                                 full.names = TRUE)




#### Read in and normalize depth data ####
raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/2018_analysis_assembly/hgcA/depth/", "", .) %>%
                               gsub("_hgcA_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        col.names = c("seqID", "depth", "length")) %>%
                               mutate(depth = depth * normalized.coverage.vector[metagenomeID]) %>%
                               mutate(read.origin = metagenomeID)
                             
                             })




#### Combine data into dataframe ####

mean.df <- do.call(rbind,
                   raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA)




#### Write out file

write.csv(mean.df,
          "dataEdited/2018_analysis_assembly/hgcA/depth/hgcA_coverage.csv",
          row.names = FALSE)
