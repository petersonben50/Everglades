#### code/2019_analysis_assembly/metabolicProteins/metabolic_proteins_depth_aggregate.R ####
# Benjamin D. Peterson


#### Clean up crew on line 5 ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(tidyverse)





#### Read in normalization vector ####

normalized.coverage.vector <- readRDS("dataEdited/2019_analysis_assembly/SCGs/scg_normalization_vector.rds")



#### Read in depth data ####

list.o.depths <- list.files(path = "dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/depth",
                            pattern = "depth.tsv",
                            full.names = TRUE)

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/depth/", "", .) %>%
                               gsub("_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        header = FALSE,
                                        col.names = c("scaffoldID", "depth")) %>%
                               mutate(depth = round(depth * normalized.coverage.vector[metagenomeID], 4)) %>%
                               mutate(read.origin = metagenomeID)
                             
                           })



#### Combine data into dataframe ####

mean.df <- do.call(rbind,
                   raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA)


colSums(mean.df[, -1])

#### Write out file

write.csv(mean.df,
          "dataEdited/2019_analysis_assembly/metabolicProteins/hydrogenases/depth/hydrogenase_depth_clean.csv",
          row.names = FALSE)
