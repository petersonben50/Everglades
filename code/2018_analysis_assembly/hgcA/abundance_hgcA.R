#### code/2018_analysis_assembly/hgcA/abundance_hgcA.R ####
# Benjamin D. Peterson

# Scripts to look at abundance of hgcA in a myriad of ways
# Can't believe I wrote "a myriad of ways"

#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(gridExtra)
library(readxl)
library(tidyverse)





#### Read in abundance and phylogenetic info ####

all.data <- full_join(read.csv("dataEdited/2018_analysis_assembly/hgcA/depth/hgcA_coverage_final.csv",
                               stringsAsFactors = FALSE),
                      read_xlsx("dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_phylogenetic_clusters.xlsx")) %>%
  arrange(assignment)



#### Set up color vector ####

CB.color.vector <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
color.vector <- c(CB.color.vector["black"],
                  CB.color.vector["orange"],
                  CB.color.vector["reddishpurple"],
                  CB.color.vector["skyblue"])
names(color.vector) <- all.data$assignment %>%
  unique() %>%
  sort()



#### Read in metadata ####

metadata <- read_xlsx("metadata/metagenomes/2018_MGs.xlsx")
metadata.vector <- metadata$siteID
names(metadata.vector) <- metadata$metagenomeID

metadata.vector[length(metadata.vector) + 1] <- "total"
names(metadata.vector)[length(metadata.vector)] <- "total"


#### Generate vector of correct order of samples along sulfate gradient ####

MG.order <- c("WCA-2A-P", "WCA-2A-A", "WCA-3A-O", "WCA-3A-F", "WW")




#### Check out overall abundances of these seqs ####

depth.of.clusters <- all.data %>%
  gather(key = siteID,
         value = coverage,
         c(4:9)) %>%
  mutate(siteID = metadata.vector[siteID]) %>%
  group_by(siteID, assignment) %>%
  summarize(coverage = round(sum(coverage), 3),
            count = n()) %>%
  spread(key = siteID,
         value = coverage) %>%
  select(assignment, count, all_of(MG.order), total) 





#### First focus on summary stats ####

summary.stats.hgcA <- depth.of.clusters %>%
  select(assignment, count, total) %>%
  mutate(`Percent of total\ncoverage` = round(total / sum(depth.of.clusters$total) * 100, 2)) %>%
  rename(`Taxonomic Group` = assignment) %>%
  rename(`Number of\nSequences` = count)

# Write out a nice table with that data
png("results/2018_analysis_assembly/hgcA/depth/taxonomic_cluster_summary_stats.png",
    width = 5,
    height = 2,
    unit = "in",
    res = 150)
grid.table(summary.stats.hgcA,
           rows = NULL)
dev.off()




#### Look at depth at each site ####

coverage.per.site <- depth.of.clusters %>%
  select(assignment, all_of(MG.order)) %>%
  rename(`Taxonomic Group` = assignment)


# Relative depth
rel.coverage.per.site <- coverage.per.site
rel.coverage.per.site[, -1] <- apply(rel.coverage.per.site[, -1],
                                     1,
                                     function(x) {
                                       x / colSums(rel.coverage.per.site[, -1]) * 100
                                     }) %>% t() %>% as.data.frame()





# Plot overall and relative abundance together

png("results/2018_analysis_assembly/hgcA/depth/depth_by_taxonomic_cluster_plot.png",
    width = 8,
    height = 8,
    unit = "in",
    res = 150)

par(mar = c(5, 3, 1, 1),
    mfrow = c(2, 1),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
barplot(coverage.per.site[, -1] %>% as.matrix(),
        col = color.vector[coverage.per.site$`Taxonomic Group`],
        xlim = c(0, 6),
        ylim = c(0, 14),
        space = 0.1,
        width = 0.9,
        xaxt = "n",
        yaxt = "n",
        border = "black",
        las = 2)
axis(1,
     at = c(1:5) - 0.5,
     labels = names(coverage.per.site)[-1],
     tick = FALSE,
     las = 2,
     cex.axis = 0.8)
axis(2, 
     at = seq(0, 14, by = 2),
     las = 2,
     tck = -0.01)
title(ylab = "hgcA coverage",
      line = 2)
legend(x = 5,
       y = 10,
       legend = names(color.vector),
       fill = color.vector,
       cex = 0.8)

barplot(rel.coverage.per.site[, -1] %>% as.matrix(),
        col = color.vector[rel.coverage.per.site$`Taxonomic Group`],
        xlim = c(0, 6),
        ylim = c(0, 100),
        space = 0.1,
        width = 0.9,
        xaxt = "n",
        yaxt = "n",
        border = "black",
        las = 2)
axis(1,
     at = c(1:5) - 0.5,
     labels = names(rel.coverage.per.site)[-1],
     tick = FALSE,
     las = 2,
     cex.axis = 0.8)
axis(2, 
     at = seq(0, 100, by = 20),
     las = 2,
     tck = -0.01)
title(ylab = "Relative hgcA coverage",
      line = 2)

dev.off()



#### Look at abundance patterns of individual sequences ####

plot.individual.seq.coverage <- function(sequence.of.interest) {
  
  all.data.seq <- all.data %>%
    filter(seqID == sequence.of.interest) %>%
    select(-c(scaffoldID, seqID, length, total))
  
  taxonomic.group <- all.data.seq$assignment
  
  barplot(all.data.seq[, 1:5] %>% as.matrix(),
          col = color.vector[taxonomic.group],
          xlim = c(0, 5),
          ylim = c(0, 2.7),
          xaxt = "n",
          xaxs = "i",
          space = 0.1,
          width = 0.9,
          border = "black",
          las = 2)
  
  title(ylab = "hgcA coverage",
        line = 2.25)
  title(main = sequence.of.interest,
        cex.main = 0.8)
  
}



png("results/2018_analysis_assembly/hgcA/depth/individual_seq_coverage.png",
    width = 9,
    height = 9,
    unit = "in",
    res = 150)
par(mar = c(1, 3.5, 1, 1),
    mfrow = c(5, 5),
    mgp=c(1.5,0.4,0),
    tck=-0.008)
sapply(all.data$seqID,
       plot.individual.seq.coverage)
dev.off()