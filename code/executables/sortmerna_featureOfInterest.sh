#!/bin/sh


##############################
# code/metagenomeBinProcessing/anvioDB_generation.sh
# Benjamin D. Peterson

# This is a general script to read in
# a list of bins, a folder of metagenomes,
# and generate individual anvioDBs for each
# bin.
##############################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate sortmerna
PYTHONPATH=""

cd $workingDirectory
echo "Mapping against these databases:"
ls -l $reference

sortmerna --ref $reference \
          --reads $MG_subset \
          --fastx \
          -v \
          --paired_in \
          --threads 3 \
          --workdir sortmernaWD_$MG_subset/ \
          --aligned $featureName\_mappingTo \
          --other $featureName\_nonMapping
rm -rf sortmernaWD_$MG_subset
