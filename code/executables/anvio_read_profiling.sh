#!/bin/sh

#########################
# anvio_read_profiling.sh
# Benjamin D. Peterson

# This script takes the reads from each
# metagenome that have been mapped to the
# assemblies from this analysis.
#

# The submission file loops over each assembly.
#########################


# Will read in these variables from
# submission file:
# $assembly

# Will read in these paths from submission file:
# metagenomeList: List of metagenomes that were mapped to the assembly for binning
# output: Output location. Should be same place as the anvio databases
# mapping: Location of mapping files


#########################
# Activate environment
#########################
cd $output
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""


#########################
# Generate read profiles
#########################
if [ ! -d $assembly.merged ]; then
  cat $metagenomeList | while read metagenome
  do
    if [ ! -d $metagenome\_to_$assembly.profile ]; then
      echo "Profiling reads mapped from:" $metagenome "to" $assembly
      anvi-profile -c $assembly.db \
                    -i $mapping/$metagenome\_to_$assembly.bam \
                    --write-buffer-size 500 \
                    --min-contig-length 2000 \
                    --num-threads 10 \
                    -o $metagenome\_to_$assembly.profile
    else
      echo "STOP: Profiling reads mapped to:" $assembly "from" $metagenome "is already done"
    fi
  done
  anvi-merge *to_$assembly.profile/PROFILE.db \
              -o $assembly.merged \
              -c $assembly.db \
              -S $assembly\_merged \
              --skip-hierarchical-clustering
else
  echo "Merging profiles complete for" $groupName
fi
