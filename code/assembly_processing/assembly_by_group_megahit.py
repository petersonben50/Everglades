#################################
# coassembly_by_group_megahit.py
# Benjamin D. Peterson

# This script will read in a csv file with the
# group information for coassembling metagenomes
# and start megahit on these groups.
#################################

#######################################
# Load needed python packages
#######################################
import os
import sys
import numpy
import pandas as pd

#######################################
# Read in the needed inputs from the command line
#######################################
group = sys.argv[1]
clusterSpreadSheetName = sys.argv[2]
readLocation = sys.argv[3]
output = sys.argv[4]

#######################################
# Read in the grouping file
#######################################
group2metagenome = pd.read_csv(clusterSpreadSheetName)

#######################################
# Select the metagenomes that correspond to the selected group
#######################################
groupMetagenomes = group2metagenome.loc[group2metagenome['assemblyID'] == group]
neededMetagenomes = groupMetagenomes['metagenomeID']

#######################################
# Loop over the metagenome names, adding them to the MegaHit command
#######################################
megahitCommand = 'megahit --k-min 21 --k-max 121 --k-step 10 -t 16 '
for metagenome in neededMetagenomes:
    print(metagenome)
    megahitCommand = megahitCommand + '-1 ' + readLocation + '/' + metagenome + '_R1.fastq.gz '
    megahitCommand = megahitCommand + '-2 ' + readLocation + '/' + metagenome + '_R2.fastq.gz '
    megahitCommand = megahitCommand + '-r ' + readLocation + '/' + metagenome + '_single.fastq.gz '
    print(megahitCommand)

megahitCommand = megahitCommand + '-o ' + output

#######################################
# Run megahit
#######################################
os.system(megahitCommand)
