rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(readxl)
library(tidyverse)



##### Read in MG metadata ####
MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.renaming.vector <- MG.metadata$siteID
names(site.renaming.vector) <- MG.metadata$metagenomeID

medium.renaming.vector <- MG.metadata$medium
names(medium.renaming.vector) <- MG.metadata$metagenomeID


#### Read in data ####
hgcA.read.counts <- read.table("dataEdited/2019_readBasedAnalysis/hgcA/hgcA_reads_count.txt",
                               header = TRUE) %>%
  mutate(hgcA_percentage = hgcA_reads / (non_hgcA_reads + hgcA_reads) * 100) %>%
  mutate(siteID = site.renaming.vector[metagenomeID]) %>%
  mutate(medium = medium.renaming.vector[metagenomeID])


#### Set site order for plotting ####
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")



#### Make plot ####
hgcA.read.counts %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  filter(medium == "sediment")%>%
  ggplot(aes(x = siteID,
             y = hgcA_percentage)) +
  geom_point()

# There is clearly an error in how sortmerna is processing the reads,
# as it's a bit absurd that hgcA accounts for nearly 0.1% of the reads
# in some samples, as I outlined in my Obsidian notebook on 2021-02-02.

# Let's leave this analysis alone for now. I think this approach could 
# be valuable and may merit additional analysis in the future, but not
# right now. There may be some errors in the fasta files. Doing some
# mapping of the so-called mapped reads to the database may prove
# enlightening.