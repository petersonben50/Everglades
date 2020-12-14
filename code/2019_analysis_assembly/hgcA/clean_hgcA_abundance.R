#### code/2019_analysis_assembly/clean_hgcA_coverage.R ####
# Written by Benjamin D. Peterson


#### First: Normalize to the read count ####

#### Get ready. Get set ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(stringr)
library(tidyverse)



#### Read in normalization vector ####

normalized.coverage.vector <- readRDS("dataEdited/metagenomes/reports/metagenome_normalization_vector.rds")



# Read in list of good hgcA seqs
hgcA.list <- readLines("dataEdited/2019_analysis_assembly/hgcA/identification/hgcA_good.txt")



#### List of file names with depth info ####

list.o.depth.files <- list.files(path = "dataEdited/2019_analysis_assembly/hgcA/depth",
                                 pattern = "hgcA_depth",
                                 full.names = TRUE)




#### Read in and normalize depth data ####
raw.depth.counts <- lapply(list.o.depth.files,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/2019_analysis_assembly/hgcA/depth/", "", .) %>%
                               gsub("_hgcA_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        col.names = c("seqID", "depth", "length")) %>%
                               mutate(depth = depth * normalized.coverage.vector[metagenomeID]) %>%
                               mutate(read.origin = metagenomeID)
                             
                             })





#### Combine data into dataframe ####

depth.data <- do.call(rbind,
                      raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA) %>%
  rename(scaffoldID = seqID)






#### Write out file with depths for all sequences ####

write.csv(depth.data,
          "dataEdited/2019_analysis_assembly/hgcA/depth/hgcA_coverage_all.csv",
          row.names = FALSE)






#### Clean depth file for final hgcA seqs ####

# Get a scaffold to ORF vector
scaffold.to.ORF <- sapply(depth.data$scaffoldID,
                          function(x) {
                            grep(x, hgcA.list, value = TRUE)
                          }) %>%
  unlist()

# Keep only depth info for good seqs and sum abundance
depth.data.final <- depth.data %>%
  mutate(seqID = scaffold.to.ORF[scaffoldID]) %>%
  filter(!is.na(seqID)) %>%
  mutate(total = KMBP005A + KMBP005B + KMBP005C + KMBP005D + KMBP005E +
           KMBP005F + KMBP005G + KMBP005H + KMBP005I + KMBP005J + KMBP005K +
           KMBP005L + KMBP006A + KMBP006B + KMBP006C + KMBP006D + KMBP006E +
           KMBP006F) %>%
  select(seqID, scaffoldID, length, total, KMBP005A, KMBP005B,
         KMBP005C, KMBP005D, KMBP005E, KMBP005F, KMBP005G,
         KMBP005H, KMBP005I, KMBP005J, KMBP005K, KMBP005L,
         KMBP006A, KMBP006B, KMBP006C, KMBP006D, KMBP006E, 
         KMBP006F)





#### Read out depth file for final hgcA seqs ####
write.csv(depth.data.final,
          "dataEdited/2019_analysis_assembly/hgcA/depth/hgcA_coverage.csv",
          row.names = FALSE)

















#### Second: Normalize to the SCG coverage ####

#### Get ready. Get set ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(stringr)
library(tidyverse)



#### Read in normalization vector ####

normalized.coverage.vector <- readRDS("dataEdited/2019_analysis_assembly/SCGs/scg_normalization_vector.rds")



# Read in list of good hgcA seqs
hgcA.list <- readLines("dataEdited/2019_analysis_assembly/hgcA/identification/hgcA_good.txt")



#### List of file names with depth info ####

list.o.depth.files <- list.files(path = "dataEdited/2019_analysis_assembly/hgcA/depth",
                                 pattern = "hgcA_depth",
                                 full.names = TRUE)




#### Read in and normalize depth data ####
raw.depth.counts <- lapply(list.o.depth.files,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/2019_analysis_assembly/hgcA/depth/", "", .) %>%
                               gsub("_hgcA_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        col.names = c("seqID", "depth", "length")) %>%
                               mutate(depth = depth * normalized.coverage.vector[metagenomeID]) %>%
                               mutate(read.origin = metagenomeID)
                             
                           })





#### Combine data into dataframe ####

depth.data <- do.call(rbind,
                      raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA) %>%
  rename(scaffoldID = seqID)






#### Write out file with depths for all sequences ####

write.csv(depth.data,
          "dataEdited/2019_analysis_assembly/hgcA/depth/hgcA_coverage_all_scgNormalization.csv",
          row.names = FALSE)






#### Clean depth file for final hgcA seqs ####

# Get a scaffold to ORF vector
scaffold.to.ORF <- sapply(depth.data$scaffoldID,
                          function(x) {
                            grep(x, hgcA.list, value = TRUE)
                          }) %>%
  unlist()

# Keep only depth info for good seqs and sum abundance
depth.data.final <- depth.data %>%
  mutate(seqID = scaffold.to.ORF[scaffoldID]) %>%
  filter(!is.na(seqID)) %>%
  mutate(total = KMBP005A + KMBP005B + KMBP005C + KMBP005D + KMBP005E +
           KMBP005F + KMBP005G + KMBP005H + KMBP005I + KMBP005J + KMBP005K +
           KMBP005L + KMBP006A + KMBP006B + KMBP006C + KMBP006D + KMBP006E +
           KMBP006F) %>%
  select(seqID, scaffoldID, length, total, KMBP005A, KMBP005B,
         KMBP005C, KMBP005D, KMBP005E, KMBP005F, KMBP005G,
         KMBP005H, KMBP005I, KMBP005J, KMBP005K, KMBP005L,
         KMBP006A, KMBP006B, KMBP006C, KMBP006D, KMBP006E, 
         KMBP006F)





#### Read out depth file for final hgcA seqs ####
write.csv(depth.data.final,
          "dataEdited/2019_analysis_assembly/hgcA/depth/hgcA_coverage_scgNormalization.csv",
          row.names = FALSE)

