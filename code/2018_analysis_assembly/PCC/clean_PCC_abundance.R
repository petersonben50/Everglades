#### code/2018_analysis_assembly/PCC/clean_PCC_abundance ####
# Benjamin D. Peterson


#### Clean up crew on line 7 ####

rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(tidyverse)



#### Read in normalization vector ####

normalized.coverage.vector <- readRDS("dataEdited/metagenomes/reports/metagenome_normalization_vector.rds")



#### Read in depth data ####

list.o.depths <- list.files(path = "dataEdited/2018_analysis_assembly/metabolicProteins/PCC/BBOMP_depth",
                            pattern = "depth.tsv",
                            full.names = TRUE)

raw.depth.counts <- lapply(list.o.depths,
                           function(fileName) {
                             
                             metagenomeID = fileName %>%
                               gsub("dataEdited/2018_analysis_assembly/metabolicProteins/PCC/BBOMP_depth/", "", .) %>%
                               gsub("_depth.tsv", "", .)
                             
                             read.table(fileName,
                                        stringsAsFactors = FALSE,
                                        header = FALSE,
                                        col.names = c("scaffoldID", "depth")) %>%
                               mutate(depth = depth * normalized.coverage.vector[metagenomeID]) %>%
                               mutate(read.origin = metagenomeID)
                             
                           })



#### Combine data into dataframe ####

mean.df <- do.call(rbind,
                   raw.depth.counts) %>%
  spread(key = read.origin, value = depth, fill = NA)




#### Write out file ####

write.csv(mean.df,
          "dataEdited/2018_analysis_assembly/metabolicProteins/depth/PCC_depth_clean.csv",
          row.names = FALSE)




#### Generate plot ####
# Read in metadata
metadata <- read_xlsx("metadata/metagenomes/2018_MGs.xlsx")
metadata.vector <- metadata$siteID
names(metadata.vector) <- metadata$metagenomeID

# Generate vector of correct order of samples along sulfate gradient
MG.order <- c("WCA-2A-P", "WCA-2A-A", "WCA-3A-O", "WCA-3A-F", "WW")

# Plotting function
plot.scaffold.coverage <- function(scaffold.list,
                                   geneName) {
  
  mean.df %>%
    filter(scaffoldID %in% scaffold.list) %>%
    gather(key = sampleID,
           value = coverage,
           -1) %>%
    mutate(sampleID = metadata.vector[sampleID]) %>%
    mutate(sampleID = fct_relevel(sampleID, MG.order)) %>%
    ggplot(aes(y=coverage, x=sampleID)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = geneName) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
}



plot.scaffold.coverage(mean.df$scaffoldID,
                       "PCC")
