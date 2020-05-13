#!/bin/sh

######################
# code/2018_analysis_assembly/sulfur/sulfur_genes.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Sulfate reducers
############################################
############################################

######################
# dsrA
######################

screen -S EG_dsrA_tree
scripts=~/Everglades/code/generalUse
# Set up directory
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins
mkdir sulfur
mkdir sulfur/dsr
mkdir sulfur/dsr/dsrA_phylogeny

# Copy over my sequences
cp derep/dsrA.afa sulfur/dsr/dsrA_phylogeny/dsrA.afa
# Copy in reference sequences
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny/
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny/
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa


# Generate ML tree
screen -S EG_dsrA_tree
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_phylogeny/
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s dsrA_phylogeny_masked.afa \
        -n dsrA
