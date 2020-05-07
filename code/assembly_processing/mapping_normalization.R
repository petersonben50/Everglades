#### code/assembly_processing/mapping_normalization.R ####
# Benjamin D. Peterson

# This script will read in the coverage of the
# metagenomes from the Everglades. It will read
# out a vector for normalizing the coverages of
# sequences in those metagenomes.



#### Get your junk out of here ####

rm(list = ls())
setwd("Documents/research/Everglades/")
library(tidyverse)



#### Read in data ####

coverage.data <- read.table("dataEdited/metagenomes/reports/metagenome_coverage.tsv",
                            sep = '\t',
                            header = TRUE,
                            stringsAsFactors = FALSE)









#### Normalization value ####

# Set the values used in normalizing coverage

# Assumed length of reads:
read_length <- 150

# Paired end or not?
# Paired = 2
# Single = 1
number_of_ends <- 2

# Number of reads to use:
# (We'll normalize to 10 million reads).
read_count <- 10000000

# Final coverage for normalization
normal_coverage <- read_length * number_of_ends * read_count

# Clean up
rm(read_length, number_of_ends, read_count)










#### Calculate normalization values ####

# Sum the coverage from the different types
# of reads
coverage.data <- coverage.data %>%
  mutate(total_coverage = R1 + R2 + single + merged)

# Generate vector
coverage.vector <- coverage.data$total_coverage
names(coverage.vector) <- coverage.data$metagenomeID

# Normalize vector
normalized.coverage.vector <- normal_coverage / coverage.vector

# Clean up
rm(coverage.data, coverage.vector, normal_coverage)


#### Write out R object of the normalized vector ####

saveRDS(object = normalized.coverage.vector,
        file = "dataEdited/metagenomes/reports/metagenome_normalization_vector.rds")
