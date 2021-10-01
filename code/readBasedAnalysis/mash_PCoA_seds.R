#### code/mash/NMDS_plot_MGs_2019_seds.R ####
# Benjamin D. Peterson

#### Always start with a clean slate ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(ape)
library(lubridate)
library(patchwork)
library(readxl)
library(rgl)
library(scatterplot3d)
library(tidyverse)
library(vegan)
library(vegan3d)
cb.translator <- readRDS("/Users/benjaminpeterson/Box/ancillary_science_stuff/colors/colorblind_friendly_colors_R/colorblind_friendly_colors.rds")
MG.order <- c("2A-N", "2A-A", "3A-O", "3A-N", "3A-F", "LOX8")



#### Read in distance data ####

dist.data <- read.table("dataEdited/metagenomes/mash_data/EG_MG_2019_seds.dist",
                        stringsAsFactors = FALSE,
                        sep = "\t")
names(dist.data) <- c("refID", "queryID", "mash_dist",
                      "pvalue", "matching_hashes")
clean.dist.data <- dist.data %>%
  mutate(refID = refID %>%
           gsub("temp_MG_files/", "", .) %>%
           gsub(".fastq.gz", "", .)) %>%
  mutate(queryID = queryID %>%
           gsub("temp_MG_files/", "", .) %>%
           gsub(".fastq.gz", "", .)) %>%
  select(refID, queryID, mash_dist) %>%
  spread(key = queryID,
         value = mash_dist)



#### Read in metadata ####

MG.data <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")
site.vector <- MG.data$siteID
names(site.vector) <- MG.data$metagenomeID



#### Generate matrix ####

clean.dist.data.matrix <- clean.dist.data %>%
  select(-refID)
row.names(clean.dist.data.matrix) <- clean.dist.data$refID



#### Calculate dimensions ####

pcoa.analysis <- pcoa(clean.dist.data.matrix)
pcoa.analysis
summary(pcoa.analysis)



# Transfer coordinates to dataframe
coordinates.df <- as.data.frame(list(pcoa.analysis$vectors[, "Axis.1"],
                                     pcoa.analysis$vectors[, "Axis.2"],
                                     pcoa.analysis$vectors[, "Axis.3"]),
                               stringsAsFactors = FALSE)
names(coordinates.df) <- c("pcoa_1", "pcoa_2", "pcoa_3")



#### Add metadata ####

coordinates.df$site <- site.vector[rownames(coordinates.df)]



#### Determine how much variance each axis accounts for ####

PCoA1_var <- round(pcoa.analysis$values[["Eigenvalues"]][1]/sum(pcoa.analysis$values[["Eigenvalues"]])*100, 1)
PCoA2_var <- round(pcoa.analysis$values[["Eigenvalues"]][2]/sum(pcoa.analysis$values[["Eigenvalues"]])*100, 1)
PCoA3_var <- round(pcoa.analysis$values[["Eigenvalues"]][3]/sum(pcoa.analysis$values[["Eigenvalues"]])*100, 1)



#### Make plots ####

# Plot axis 1 vs. axis 2
AX1_AX2 <- ggplot(coordinates.df, aes(x = pcoa_1, y = pcoa_2)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = paste("PCoA1 ", PCoA1_var, "% of variance",
                 sep = ""),
       y = paste("PCoA2 ", PCoA2_var, "% of variance",
                 sep = "")) +
  geom_text(aes(label = site)) +
  labs(title = "PCoA ordination of 2019 MGs by\nMash distance: axes 1 and 2") 

# Plot axis 2 vs. axis 3
AX2_AX3 <- ggplot(coordinates.df, aes(x = pcoa_2, y = pcoa_3)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = paste("PCoA2 ", PCoA2_var, "% of variance",
                 sep = ""),
       y = paste("PCoA3 ", PCoA3_var, "% of variance",
                 sep = "")) +
  geom_text(aes(label = site)) +
  labs(title = "PCoA ordination of 2019 MGs by\nMash distance: axes 2 and 3") 
  

#### Save out plots ####
pdf("results/metagenomes/mash/PCoA_plot_MGs_2019_seds_AX1AX2AX3.pdf",
    width = 10,
    height = 5)
AX1_AX2 + AX2_AX3
dev.off()

# Lox is really separated out from the other microbial communities, primarily on the first axis.
# This axis accounts for 32.8% of the variance.
# 3A-F and 2A-N are both divergent from the other WCA sites, and do not cluster together (in fact
# they're divergent from each other).
# 20.8% of the variance here.
# 3A-O, 2A-A, and 2A-N cluster tightly by the first two axes. However, on the third axis of variance,
# 2A-A separates out substantially from 3A-O and 3A-N



#### Save out R object to put into Rmd ####
# saveRDS(object = pcoa.analysis,
#         file = "dataEdited/metagenomes/mash_data/PCoA_2019_sed.rds")


#### Plot 3D image ####

color.vector <- cb.translator[1:6]
names(color.vector) <- MG.order



plot3d(pcoa.analysis$vectors[, 1:3],
       xlab = paste("PCoA1 ", PCoA1_var, "% of variance",
                    sep = ""),
       ylab = paste("PCoA2 ", PCoA2_var, "% of variance",
                    sep = ""),
       zlab = paste("PCoA3 ", PCoA3_var, "% of variance",
                    sep = ""),
       col = color.vector[site.vector[row.names(pcoa.analysis$vectors)]],
       size = 10)
decorate3d(main = "First 3 axes of PCoA analysis of 2019 MGs")
legend3d("topright",
         legend = names(color.vector),
         text.col = color.vector)



#### Plot 3D ordination using scatterplot3d ####
pdf("results/metagenomes/mash/PCoA_plot_MGs_2019_seds_3D.pdf",
    height = 6,
    width = 6)
par(mfrow = c(1, 1),
    mar = c(3, 3, 1, 1))
scatterplot3d(-pcoa.analysis$vectors[, 1:3],
              pch = 16,
              type = "h",
              color = color.vector[site.vector[row.names(pcoa.analysis$vectors)]],
              lab = c(4, 4),
              lab.z = 4,
              xlab = paste("PCoA1 ", PCoA1_var, "% of variance",
                           sep = ""),
              ylab = paste("PCoA2 ", PCoA2_var, "% of variance",
                           sep = ""),
              zlab = paste("PCoA3 ", PCoA3_var, "% of variance",
                           sep = ""))
legend(x = 5,
       y = 4,
       legend = names(color.vector),
       text.col = color.vector,
       bg = "white")
dev.off()




#### Run PERMANOVA on dataset ####
MG.data.permanova <- MG.data %>%
  filter(grepl("KMBP005",
               metagenomeID))
adonis2(clean.dist.data.matrix ~ sampleID,
        data = MG.data.permanova)