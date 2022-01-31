#### code/binning/aggregate_binning_information.R ####
# Benjamin D. Peterson


#### Clean up on aisle R ####
rm(list = ls())
setwd("/Users/benjaminpeterson/Documents/research/Everglades")
library(readxl)
library(tidyverse)


#### Read in hgcA lists ####
hgcA.info <- read_xlsx("dataEdited/bin_analysis/hgcA/notes_on_bin_hgcA_seqs.xlsx")


#### Add in taxonomy data ####
tax.data <- read.table("dataEdited/bin_analysis/binning_initial/binsGood/taxonomy_summary.txt",
                       sep = '\t',
                       col.names = c("binID", "taxonomy"))
hgcA.info <- hgcA.info %>%
  left_join(tax.data)
rm(tax.data)


#### Read in completeness data from anvio ####
anvio.stats <- read.table("dataEdited/bin_analysis/binning_initial/binsRaw/hgcA_bins_summary.txt",
                          sep = "\t",
                          header = TRUE)
names(anvio.stats) <- c("binID", "taxon", "genomeLength_bp_anvio", "scaffolds_anvio",
                        "N50_anvio", "GC_anvio", "Com_anvio", "Red_anvio")
hgcA.info <- hgcA.info %>%
  left_join(anvio.stats %>%
              select(binID, Com_anvio, Red_anvio))
rm(anvio.stats)


#### Read in completeness data from checkM ####
checkm.stats <- read.csv("dataEdited/bin_analysis/binning_initial/binsRaw/checkM_stats.csv")
names(checkm.stats) <- c("binID", "Com_checkm", "Red_checkm", "strainHet",
                         "genomeSize_bp", "scaffolds", "N50", "scaffoldMeanLength_bp",
                         "longestScaffold_bp", "GC", "predictedGenes")
hgcA.info <- hgcA.info %>%
  left_join(checkm.stats)
rm(checkm.stats)


#### Read in coverage data ####

cov.data <- read.table("dataEdited/bin_analysis/binning_initial/binsGood/coverage_goodBins.txt",
                       sep = '\t',
                       header = TRUE)
names(cov.data)[1] <- "binID"
names(cov.data)[-1] <- names(cov.data)[-1] %>%
  strsplit("_TO_") %>%
  sapply("[", 1)
hgcA.info <- hgcA.info %>%
  left_join(cov.data)
rm(cov.data)


#### Read in hgcA representative key ####
hgcA.key <- readRDS("dataEdited/assembly_analysis/hgcA/dereplication/hgcA_derep_key.rds")
hgcA.info <- hgcA.info %>%
  left_join(hgcA.key)
rm(hgcA.key)


#### Save out data ####
write.csv(hgcA.info,
          file = "dataEdited/bin_analysis/hgcA_bin_data.csv",
          row.names = FALSE)
