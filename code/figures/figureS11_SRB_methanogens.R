#### code/figures/figureS9_SRB_methanogens.R ####
# Benjamin D. Peterson


#### Rock em sock em ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades/")
library(ggpubr)
library(tidyverse)


#### Read in ggplots of needed graphs ####
dsr.abund <- readRDS("results/metagenomes/assembly/SRB/dsr_abundance.rds")
dsr.ord <- readRDS("results/metagenomes/assembly/SRB/dsrA_ordination.rds")
mcrA.abund <- readRDS("results/metagenomes/assembly/methanogenesis/mcrA_abundance.rds")
mcrA.ord <- readRDS("results/metagenomes/assembly/methanogenesis/mcrA_ordination.rds")


#### Save out figure ####
pdf("results/figures/S9_SRB_methanogens.pdf",
    width = 10,
    height = 7.5)
ggarrange(dsr.abund, dsr.ord,
          mcrA.abund, mcrA.ord,
          ncol = 2,
          nrow = 2,
          widths = c(1, 1),
          labels = c("A.", "B.", "C.", "D."))
dev.off()

