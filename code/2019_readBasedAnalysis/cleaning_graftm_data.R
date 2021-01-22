#### Clean up first ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(tidyverse)


#### Read in files ####
list.o.files <- list.files("dataEdited/2019_readBasedAnalysis/16S/graftM/count_table",
                           pattern = ".txt")
files.of.interest <- list(0)
for (file.number in 1:length(list.o.files)) {
  files.of.interest[[file.number]] <- read.table(paste("dataEdited/2019_readBasedAnalysis/16S/graftM/count_table/",
                                                       list.o.files[file.number],
                                                       sep = ""),
                                                 sep = '\t',
                                                 header = TRUE,
                                                 comment.char="") %>%
    select(-X.ID)
}



#### Join all dataframes ####
taxonomy.df <- Reduce(function(x, y, ...) { full_join(x, y) },
                      files.of.interest) %>%
  gather(value = counts,
         key = metagenomeID,
         -2) %>%
  mutate(metagenomeID = metagenomeID %>% gsub("_R1", "", .))



#### Relativize data ####
summed.counts <- taxonomy.df %>%
  group_by(metagenomeID) %>%
  summarise(total.counts = sum(counts, na.rm = TRUE))
taxonomy.df <- taxonomy.df %>%
  left_join(summed.counts) %>%
  mutate(rel.abundance = counts / total.counts) %>%
  select(metagenomeID, ConsensusLineage, counts, rel.abundance)

#### Replace eukaryotes names with just kingdom ####
taxonomy.df[grep("Root; Eukaryota", x = taxonomy.df$ConsensusLineage), "ConsensusLineage"] <- "Root; Eukaryota"



#### Split out taxonomy ####
taxonomy.df.split <- taxonomy.df %>%
  separate(col = ConsensusLineage,
           into = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
           "; ") %>%
  select(-root)


#### Save out data ####
saveRDS(taxonomy.df.split,
        "dataEdited/2019_readBasedAnalysis/16S/graftM/taxonomy_counts.rds")
