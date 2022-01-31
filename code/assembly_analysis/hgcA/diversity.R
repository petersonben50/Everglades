#### code/assembly_analysis/hgcA/diversity.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggpubr)
library(phyloseq)
library(readxl)
library(scatterplot3d)
library(tidyverse)
library(vegan)
source("code/setup_PW_core_order_color_points.R")


#### Read in metadata ####
metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx") %>%
  filter(grepl("KMBP005",
               x = metagenomeID))
# Site metadata
site.metadata.vector <- metadata$siteID
names(site.metadata.vector) <- metadata$metagenomeID
# Medium metadata
medium.metadata.vector <- metadata$medium
names(medium.metadata.vector) <- metadata$metagenomeID


#### Set up color vector for microbes ####
color.code.df <- read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                           sheet = "color_codes")
color.code.vector.fused <- cb.translator[color.code.df$colorCode]
names(color.code.vector.fused) <- color.code.df$clusterName
color.code.vector <- color.code.vector.fused[which(names(color.code.vector.fused) != "fused")]


#### Read in hgcA abundance and phylogenetic info ####
all.data <- read.csv("dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv",
                     stringsAsFactors = FALSE) %>%
  full_join(read_xlsx("dataEdited/assembly_analysis/hgcA/phylogeny/seq_classification.xlsx",
                      sheet = "seq_class")) %>%
  arrange(classification) %>%
  filter(seqID %in% readLines("dataEdited/assembly_analysis/hgcA/hgcA_true.txt")) %>%
  filter(clusterName != "fused") %>%
  filter(clusterName != "paralog") %>%
  gather(key = MG,
         value = coverage,
         c(4:22)) %>%
  filter(MG != "total") %>%
  mutate(siteID = site.metadata.vector[MG],
         medium = medium.metadata.vector[MG]) %>%
  filter(medium == "sediment") %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  as.data.frame()


#### Calculate "richness" ####
richness.data <- all.data %>%
  filter(coverage > 0) %>%
  group_by(MG, siteID) %>%
  summarise(unique_sequences = n())
richness.plot <- richness.data %>%
  ggplot(aes(x = siteID,
             y = unique_sequences,
             color = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  ylim(c(0, 70)) +
  ylab("Richness") +
  theme_bw() +
  theme(legend.position = "none",
        # legend.position = c(0.8, 0.4),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))
saveRDS(richness.plot,
        "results/metagenomes/assembly/hgcA/richness.rds")



#### Evenness - Shannon diversity ####
total.data <- all.data %>%
  group_by(MG) %>%
  summarise(total.coverage = sum(coverage))
evenness.data <- all.data %>%
  left_join(total.data) %>%
  mutate(relative_coverage = coverage / total.coverage) %>%
  select(MG, seqID, relative_coverage) %>%
  spread(key = seqID,
         value = relative_coverage)
row.names(evenness.data) <- evenness.data$MG
evenness.data <- evenness.data %>%
  select(-MG) %>%
  as.matrix()

# Generate phyloseq object
hgcA.phyloseq <- otu_table(evenness.data,
                           taxa_are_rows = FALSE)

# Calculate Shannon diversity
diversity.vector <- diversity(hgcA.phyloseq,
                              index = "shannon")
diversity.df <- data.frame("MG" = names(diversity.vector),
                           "diversity" = diversity.vector) %>%
  mutate(siteID = site.metadata.vector[MG])
row.names(diversity.df) <- NULL

# Plot Shannon diversity
diversity.df %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  ggplot(aes(x = siteID,
             y = diversity)) +
  geom_point() +
  ylim(c(0, 4)) +
  theme_bw()


#### Shannon's equitability ####
E.data <- diversity.df %>%
  left_join(richness.data %>% select(-siteID)) %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  mutate(E = diversity / log(unique_sequences))
E.plot <- E.data %>%
  ggplot(aes(x = siteID,
             y = E,
             color = siteID)) +
  geom_point(size = 3) +
  scale_color_manual(values = color.vector,
                     name = "Site ID") +
  ylim(c(0.5, 1)) +
  ylab("Shannon's equitability") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"))

saveRDS(E.plot,
        "results/metagenomes/assembly/hgcA/evenness.rds")


#### Plot em ####
ggarrange(richness.plot,
          E.plot)


#### Statistical analysis of evenness ####
richness.model <- aov(unique_sequences ~ siteID, data = richness.data)
summary(richness.model)
TukeyHSD(richness.model, conf.level=.95) 

evenness.model <- aov(E ~ siteID, data = E.data)
summary(evenness.model)
TukeyHSD(evenness.model, conf.level=.95) 



#### Prepare data for ordination ####

#make metadata a phyloseq class object
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata[, "metagenomeID"] %>%
  unlist()
meta.phylo <- sample_data(metadata)

##make biom for phyloseq
all.phyloseq <- merge_phyloseq(hgcA.phyloseq, meta.phylo)


#### Generate ordination ####
hgcA.bray <- ordinate(all.phyloseq,
                      method = "PCoA",
                      distance = "bray")

bc.ord.hgcA <- plot_ordination(all.phyloseq,
                               hgcA.bray,
                               color = "siteID",
                               shape = "siteID") +
  scale_color_manual(values = color.vector[1:6],
                     name = "Site ID") +
  scale_shape_manual(values = point.vector[1:6],
                     name = "Site ID") +
  geom_point(size=5) +
  ylim(c(-0.5, 0.5)) +
  xlim(c(-0.55, 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour="black"),
        axis.text.y = element_text(colour="black"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = c(0.3, 0.7))

saveRDS(bc.ord.hgcA,
        "results/metagenomes/assembly/hgcA/ordination.rds")


#### Calculate variability ####
PCoA1_var <- round(hgcA.bray$values[["Eigenvalues"]][1]/sum(hgcA.bray$values[["Eigenvalues"]])*100, 1)
PCoA2_var <- round(hgcA.bray$values[["Eigenvalues"]][2]/sum(hgcA.bray$values[["Eigenvalues"]])*100, 1)
PCoA3_var <- round(hgcA.bray$values[["Eigenvalues"]][3]/sum(hgcA.bray$values[["Eigenvalues"]])*100, 1)




#### Correlate with environmental variables ####

Hg.inc.data <- readRDS("dataEdited/incubations/incubation_data_with_normalization.rds") %>%
  group_by(coreID) %>%
  summarise(RMP_core = mean(RMP_core)) %>%
  rename(siteID = coreID)
sulfide.data <- read_xlsx("dataEdited/geochem/geochem_data_2019.xlsx",
                          sheet = "PW_SW_geochem") %>%
  mutate(siteID = fct_relevel(siteID, MG.order)) %>%
  filter(medium == "PW") %>%
  select(siteID, Sulfide_µg.L)
data.matrix <- full_join(Hg.inc.data,
                         sulfide.data) %>%
  arrange(siteID) %>%
  select(-siteID) %>%
  as.matrix()
row.names(data.matrix) <- Hg.inc.data$siteID

envfit(hgcA.bray, data.matrix)

#### Generate 3D ordination ####
scatterplot3d(hgcA.bray$vectors[, 1:3],
              pch = 16,
              type = "h",
              color = color.vector[site.metadata.vector[row.names(hgcA.bray$vectors)]],
              lab = c(4, 4),
              lab.z = 4,
              # zlim = c(-0.3, 0.7),
              # ylim = c(-0.3, 0.7),
              # xlim = c(-0.6, 0.6),
              xlab = paste("PCoA1: ", PCoA1_var, "% of variance",
                           sep = ""),
              ylab = paste("PCoA2: ", PCoA2_var, "% of variance",
                           sep = ""),
              zlab = paste("PCoA3: ", PCoA3_var, "% of variance",
              sep = ""),
)
legend(x = 4.5,
       y = 5.5,
       legend = names(color.vector),
       text.col = color.vector,
       bg = "white")

