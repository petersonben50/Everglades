#### code/2019_binning/binning_stats_dereplication.R ####
# Benjamin D. Peterson

# This script will aggregate the bin information and
# save out a spreadsheet with that information that
# I can use for manual dereplication.


#### Set the table ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(igraph)
library(tidyverse)
library(vegan)

#### Read in completeness data ####
anvio.stats <- read.table("dataEdited/2019_binning/binning_initial/binsRaw/hgcA_bins_summary.txt",
                          sep = "\t",
                          header = TRUE)
names(anvio.stats) <- c("binID", "taxon", "genomeLength_bp_anvio", "scaffolds_anvio",
                        "N50_anvio", "GC_anvio", "Com_anvio", "Red_anvio")
anvio.stats <- anvio.stats %>%
  select(binID, Com_anvio, Red_anvio)

checkm.stats <- read.csv("dataEdited/2019_binning/binning_initial/binsRaw/checkM_stats.csv")
names(checkm.stats) <- c("binID", "Com_checkm", "Red_checkm", "strainHet",
                         "genomeSize_bp", "scaffolds", "N50", "scaffoldMeanLength_bp",
                         "longestScaffold_bp", "GC", "predictedGenes")

all.data <- full_join(anvio.stats, checkm.stats)



#### Read in taxonomy data ####

tax.data <- read.table("dataEdited/2019_binning/binning_initial/binsGood/taxonomy_summary.txt",
                       sep = '\t',
                       col.names = c("binID", "taxonomy"))
goodBins.data <- left_join(tax.data,
                           all.data)



#### Read in coverage data ####

cov.data <- read.table("dataEdited/2019_binning/binning_initial/binsGood/coverage_goodBins.txt",
                       sep = '\t',
                       header = TRUE)
names(cov.data)[1] <- "binID"
names(cov.data)[-1] <- names(cov.data)[-1] %>%
  strsplit("_TO_") %>%
  sapply("[", 1)
goodBins.data <- left_join(goodBins.data,
                           cov.data)



#### Compare C and R between programs ####

plot(x = all.data$Com_anvio,
     y = all.data$Com_checkm,
     ylab = "Completeness metric from CheckM",
     xlab = "Completeness metric from anvi'o",
     main = "Comparison of completeness data",
     pch = 18)
abline(h = 50, col = "red")
abline(v = 50, col = "red")
abline(a = 0, b = 1, col = "green")

plot(x = all.data$Red_anvio,
     y = all.data$Red_checkm,
     ylab = "Redundancy metric from CheckM",
     xlab = "Redundancy metric from anvi'o",
     main = "Comparison of redundancy data",
     pch = 18)
abline(h = 10, col = "red")
abline(v = 10, col = "red")
abline(a = 0, b = 1, col = "green")



#### Calculate HMSs ####

ANI.values <- read.table(file = 'dataEdited/2019_binning/binning_initial/binsGood/goodBins.all.ani.out.cleaned',
                         sep = "\t",
                         header = TRUE,
                         stringsAsFactors = FALSE) %>%
  mutate(GENOME1 = gsub(".fna", "", GENOME1),
         GENOME2 = gsub(".fna", "", GENOME2))
names(ANI.values)[3:6] <- c("ANI_1", "ANI_2", "AF_1", "AF_2")
ANI.values <- ANI.values %>%
  filter(ANI_1 >= 0.4)

plot(x = c(ANI.values$AF_1, ANI.values$AF_2),
     y = c(ANI.values$ANI_1, ANI.values$ANI_2),
     pch = 18,
     xlab = "Alignment Fraction",
     ylab = "ANI values")

# Set cut-offs
ANI.cutoff <- 94
cov.cutoff <- 0.3

# Generate edgelist
edgelist <- ANI.values %>%
  filter(ANI_1 > ANI.cutoff &
           ANI_2 > ANI.cutoff &
           AF_1 > cov.cutoff &
           AF_2 > cov.cutoff) %>%
  select(GENOME1, GENOME2)
# Generate the graph from the edgelist
adjacency.graph <- graph_from_data_frame(edgelist)
# Store the bin and HMS name in a df.
HMS.ID <- data.frame(paste("HMS.",
                           clusters(adjacency.graph)$membership,
                           sep = ""),
                     names(clusters(adjacency.graph)$membership),
                     stringsAsFactors = FALSE)
names(HMS.ID) <- c("HMS", "binID")
HMS.ID <- HMS.ID %>%
  mutate(binID = binID %>%
           strsplit(".fna") %>%
           sapply("[", 1)) %>%
  arrange(HMS)


# Check out loner bins
bin.list <- unique(c(ANI.values$GENOME1, ANI.values$GENOME2)) %>%
  strsplit(".fna") %>%
  sapply("[", 1)
lone.bins <- data.frame(bin.list[!(bin.list %in% HMS.ID$bin)],
                        bin.list[!(bin.list %in% HMS.ID$bin)],
                        stringsAsFactors = FALSE)
names(lone.bins) <- c("HMS", "binID")
# Add loner bins to list
HMS.bin.info <- do.call("rbind",
                        list(HMS.ID, lone.bins))
HMS.bin.vector <- HMS.bin.info$HMS
names(HMS.bin.vector) <- HMS.bin.info$binID


rm(ANI.cutoff, cov.cutoff, HMS.ID, lone.bins,
   edgelist, ANI.values, adjacency.graph,
   bin.list)





#### Check co-occurence ####
# Generate matrix for ordination
bin.depth.matrix <- cov.data %>%
  select(-binID) %>%
  as.matrix()
rownames(bin.depth.matrix) <- cov.data$binID
# bin.depth.matrix <- bin.depth.matrix[hgcA.bins.anvio, ]

# Bray-Curtis metaMDS
BC.nmds = metaMDS(bin.depth.matrix, distance="bray", k=2, trymax=1000)
data.scores = as.data.frame(scores(BC.nmds),stringsAsFactors = FALSE)
data.scores$hmsID <- HMS.bin.vector[rownames(data.scores)]

# Generate plot
ggplot(data.scores, aes(x = NMDS1,
                        y = NMDS2,
                        col = hmsID)) + 
  geom_point(size = 0)+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", y = "NMDS2") +
  geom_text(aes(label = hmsID),
            size = 2.5) +
  geom_line()



#### Add HMS info to all data ####
HMS.data.for.output <- full_join(HMS.bin.info,
                                 goodBins.data)
write.csv(HMS.data.for.output,
          "dataEdited/2019_binning/binning_initial/binsGood/HMS_data.csv",
          row.names = FALSE)
# Take this csv file and "Save As" an Excel sheet. Dereplicate the bins in one 
# of the sheets in the file. 