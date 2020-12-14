#### Clean up ####

rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(tidyverse)



#### Read in data ####

scg.abundance <- read.table("dataEdited/2019_analysis_assembly/SCGs/scg_coverage.tsv",
                            sep = "\t",
                            stringsAsFactors = FALSE,
                            col.names = c("geneID", "metagenomeID", "coverage")) 



#### Plot out raw data ####
scg.abundance %>%
  ggplot(aes(x = metagenomeID,
             y = coverage)) +
  geom_line(mapping = aes(group = geneID)) +
  geom_point() +
  theme_classic() +
  stat_summary(geom = "point", fun = "mean",
               col = "black", fill = "red",
               size = 3, shape = 24)


#### Generate normalization factor for each metagenomes ####

mean.scg.abundance <- scg.abundance %>%
  group_by(metagenomeID) %>%
  summarize(coverage = mean(coverage))
# Normalize to a coverage of 1000
normalized.mean.scg.abundance <- mean.scg.abundance %>%
  mutate(NF = 1000 / coverage)


#### Generate normalization vector for each metagenomes ####

normalization.vector <- normalized.mean.scg.abundance$NF
names(normalization.vector) <- normalized.mean.scg.abundance$metagenomeID
saveRDS(normalization.vector,
        "dataEdited/2019_analysis_assembly/SCGs/scg_normalization_vector.rds")
