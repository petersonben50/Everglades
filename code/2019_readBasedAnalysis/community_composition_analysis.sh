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
