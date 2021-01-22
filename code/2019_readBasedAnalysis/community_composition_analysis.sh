##################################################
##################################################
# code/2019_readBasedAnalysis/community_composition_analysis.sh
# Benjamin D. Peterson

# These are the bash scripts I used to analyse the
# microbial community composition using read-based
# metrics.
##################################################
##################################################

##################################################
##################################################
# Diversity analysis of metagenomes using Nonpareil
##################################################
##################################################

mkdir /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis
mkdir /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/nonpareil
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/nonpareil

chmod +x /home/GLBRCORG/bpeterson26/Everglades/code/executables/nonpareil3_run.sh
condor_submit /home/GLBRCORG/bpeterson26/Everglades/code/submission/nonpareil3_run.sub



##################################################
##################################################
# 16S analysis of metagenomes using GraftM
##################################################
##################################################

#########################
# Set up for GraftM run
#########################
mkdir /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/16S_analysis

#########################
# Submit scripts to run GraftM on all metagenomes
#########################
chmod +x /home/GLBRCORG/bpeterson26/Everglades/code/executables/graftm_run.sh
condor_submit /home/GLBRCORG/bpeterson26/Everglades/code/submission/graftm_run_silva.sub

#########################
# Clean up GraftM outputs
#########################

metagenomeList=/home/GLBRCORG/bpeterson26/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/16S_analysis
cat $metagenomeList | while read metagenome
do
  echo "Pulling out" $metagenome "16S counts"
  cp $metagenome\_output/combined_count_table.txt $metagenome\_count_table.txt
done
