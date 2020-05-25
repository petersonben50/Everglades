#### code/16S_from_metagenomes_2018/16S_abundance.R ####
# Benjamin D. Peterson



#### Clean up on aisle R ####

rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(tidyverse)



#### Read in 16S abundance ####
abundance.data <- read.table("/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/16S_from_metagenomes_2018/output/fractional_abundance_clustered.txt",
                             sep = '\t',
                             header = TRUE,
                             stringsAsFactors = FALSE) %>%
  select(seqID, metagenomeID, abundance)





#### Read in 16S taxonomy ####

taxonomy.data <- read.table("/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/16S_from_metagenomes_2018/output/16S_seqs_clean.taxonomy",
                            sep = '\t',
                            stringsAsFactors = FALSE) %>%
  select(-c("V9"))
tax.row.names <- c("seqID", "kingdom", "phylum", "class", "order", 
                   "lineage", "clade", "tribe")
names(taxonomy.data) <- tax.row.names





##### Read in MG metadata ####

MG.metadata <- read_xlsx("metadata/metagenomes/2018_MGs.xlsx")
renaming.vector <- MG.metadata$siteID
names(renaming.vector) <- MG.metadata$metagenomeID

##### Rename MGs with sample info ####

abundance.data.rename <- abundance.data %>%
  mutate(sampleID = renaming.vector[metagenomeID]) %>%
  select(-metagenomeID)

#### Set site order for plotting ####

MG.order <- c("WCA-2A-P", "WCA-2A-A", "WCA-3A-O", "WCA-3A-F", "WW")






#### Combine data ####

all.data <- full_join(abundance.data.rename,
                      taxonomy.data) %>%
  filter(kingdom != "Eukaryota")



#### Function to plot by taxa group ####
plot.taxonomic.levels <- function(taxonomic.level = "phylum",
                                  abundance.cutoff = 0.05) {
  taxa.data <- all.data
  taxa.data[, "taxa"] <- all.data[, taxonomic.level]
  taxa.data.clean <- taxa.data %>%
    group_by(taxa, sampleID) %>%
    summarise(abundance = sum(abundance)) %>%
    filter(abundance > abundance.cutoff) %>%
    ungroup() %>%
    mutate(sampleID = fct_relevel(str_trim(sampleID), MG.order))

  colorCount <- length(unique(taxa.data.clean$taxa))
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  set.seed(11)
  my_colors <- sample(getPalette(colorCount))
  names(my_colors) <-  levels(factor(unique(taxa.data.clean$taxa)))
  
  taxa.data.clean %>%
    ggplot(aes(y=abundance, x=sampleID, fill = taxa)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "Order", values = my_colors) +
    labs(title = paste("Abundance by ",
                       taxonomic.level,
                       sep = "")) +
    theme_classic() +
    theme(axis.text.x = element_text(colour="black"),
          axis.text.y = element_text(colour="black"))
  
  }





#### Plot out by group ####

png("results/16S_from_metagenomes_2018/taxonomic_levels_abundance.png",
    width = 12,
    height = 9,
    unit = "in",
    res = 200)
plot.taxonomic.levels("phylum", 0.05) / plot.taxonomic.levels("class", 0.05)
dev.off()


