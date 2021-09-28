#!/bin/sh

######################
# code/SCG/SCG_coverage_for_normalization.sh
# Benjamin D. Peterson
######################

######################
# Get set up
######################
cd /home/GLBRCORG/bpeterson26/references
awk -F ',' '{ print $1 }' /home/GLBRCORG/bpeterson26/references/scg_key.csv | \
  sort | uniq \
  > /home/GLBRCORG/bpeterson26/references/scg_list.txt
mkdir ~/Everglades/reports/
mkdir ~/Everglades/reports/scg
cd ~/Everglades/reports/scg
rm -rf errs logs outs
mkdir errs logs outs


######################
# Submit jobs
######################
cd ~/Everglades/dataEdited/2019_analysis_assembly/scg_abundance
#chmod +x /home/GLBRCORG/bpeterson26/Everglades/code/executables/SCG_abundance_in_assemblies.sh
condor_submit /home/glbrc.org/bpeterson26/Everglades/code/submission/SCG_abundance_in_assemblies_2019.sub


######################
# Concatenate files and clean up
######################
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_analysis_assembly/scg_abundance/
cat *_scg_coverage.tsv > scg_coverage.tsv
rm -rf working_directory_*
