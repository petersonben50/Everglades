#### code/2019_metagenome_info.R ####
# Benjamin D. Peterson


#### Clean up ####
rm(list = ls())
setwd("~/Documents/research/Everglades/")
library(lubridate)
library(readxl)
library(tidyverse)


#### Read in data ####

pre.trim.reads <- read.table("dataEdited/metagenomes/reports/metagenome_read_count_pre_trimming.tsv",
                             stringsAsFactors = FALSE,
                             header = TRUE) %>%
  rename(paired_reads_pre_trimming = forwardReads) %>%
  select(-reverseReads)

post.trim.reads <- read.table("dataEdited/metagenomes/reports/metagenome_read_count.tsv",
                              stringsAsFactors = FALSE,
                              header = TRUE) %>%
  rename(paired_reads = forwardReads) %>%
  select(-reverseReads)

unique.mapped.reads <- read.table("dataEdited/assembly_analysis//reports/uniq_mapped_reads_MG.tsv",
                                  stringsAsFactors = FALSE)
names(unique.mapped.reads) <- c("metagenomeID", "mappedReads_forward",
                                "mappedReads_reverse", "mappedReads_single")


coverage <- read.table("dataEdited/metagenomes/reports/metagenome_coverage.tsv",
                       stringsAsFactors = FALSE,
                       header = TRUE) %>%
  mutate(total_coverage = R1, R2, single, merged) %>%
  select(metagenomeID, total_coverage)

MG.metadata <- read_xlsx("metadata/metagenomes/metagenome_metadata.xlsx")


NP.coverage <- read.csv("results/metagenomes/housekeeping/nonpareil_coverage.csv",
                        stringsAsFactors = FALSE) %>%
  select(-siteID)



#### Calculate fraction of reads mapped ####
fraction.mapped.reads <- left_join(unique.mapped.reads,
                                   post.trim.reads) %>%
  mutate(total.mapped = mappedReads_forward + mappedReads_reverse + mappedReads_single,
         total.reads = (paired_reads * 2) + singleReads + mergedReads,
         fraction.mapped = total.mapped / total.reads * 100) %>%
  select(metagenomeID, fraction.mapped)



#### Combine data ####
metagenome.data <- MG.metadata %>%
  left_join(pre.trim.reads) %>%
  left_join(post.trim.reads) %>%
  left_join(fraction.mapped.reads) %>%
  left_join(coverage) %>%
  left_join(NP.coverage) %>%
  filter(year(ymd(date)) == 2019) %>%
  as.data.frame()



#### Write out data table ####
write.csv(x = metagenome.data,
          "results/metagenomes/housekeeping/metagenome_info.csv",
          row.names = FALSE)
