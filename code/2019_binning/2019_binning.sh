#!/bin/sh

##########################
# code/2019_binning/2019_binning.sh
# Benjamin D. Peterson

# This workflow will generate our set of
# hgcA+ genomes, as well as a set of
# uncurated bins.
##########################

####################################################
####################################################
# Prepare scaffolds and mapping files
####################################################
####################################################


##########################
# Filter out short scaffolds
##########################

screen -S EG_binning
mkdir ~/Everglades/dataEdited/2019_binning
cd ~/Everglades/dataEdited
mkdir 2019_binning/binning_initial
mkdir 2019_binning/binning_initial/scaffolds
export PATH=/home/GLBRCORG/bpeterson26/miniconda3/bin:$PATH
source activate anvio6.2
PYTHONPATH=""
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  if [ ! -e 2019_binning/binning_initial/scaffolds/$assembly\_filtered_assembly.fna ]; then
    echo "Processing" $assembly "scaffolds for binning"
      anvi-script-reformat-fasta assemblies/scaffolds/$assembly\_assembly.fna \
                              -o 2019_binning/binning_initial/scaffolds/$assembly\_filtered_assembly.fna \
                              -l 2000
  else
    echo $assembly "scaffolds already processed for binning"
  fi
done

#exit



##########################
# Map reads to filtered scaffolds
##########################

cd /home/GLBRCORG/bpeterson26/Everglades/reports/
rm -f outs/*_binningMapping.out \
      errs/*_binningMapping.err \
      logs/*_binningMapping.log

cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/
mkdir mapping
mkdir mapping/indices

cd /home/GLBRCORG/bpeterson26/Everglades/code/
chmod +x executables/binning_mapping.sh
condor_submit submission/binning_mapping.sub









####################################################
####################################################
# Generate anvio DBs
####################################################
####################################################



##########################
# Generate contig databases
##########################

screen -S EG_anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
IFS=$'\n'
cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir anvioDBs

for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  if [ ! -f anvioDBs/$assembly.db ]; then
    anvi-gen-contigs-database -f scaffolds/$assembly\_filtered_assembly.fna \
                              -o anvioDBs/$assembly.db \
                              -n $assembly
  else
    echo "Contig database for" $assembly "already exists."
  fi
done
exit




###############################
# Populate contigs DB with HMMs
###############################

screen -S EG_anvioDBs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
hgcAHMM=~/Everglades/references/hgcaAnvio
cd ~/Everglades/dataEdited/2019_binning/binning_initial
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  # Run default anvio HMMs
  anvi-run-hmms -c anvioDBs/$assembly.db \
                --num-threads 10
  # Run hgcA HMM
  anvi-run-hmms -c anvioDBs/$assembly.db \
                -H $hgcAHMM \
                --num-threads 10
done




###############################
# Add kaiju annotation information
###############################

screen -S EG_anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
kaiju=~/references/kaiju_db
cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir contigsTaxonomy/
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do

  if [ ! -e contigsTaxonomy/$assembly\_geneCalls_kaiju.names ]; then

    echo "Finding scaffold taxonomies for" $assembly

    conda activate anvio6.2
    PYTHONPATH=""
    anvi-get-sequences-for-gene-calls -c anvioDBs/$assembly.db \
                                      -o contigsTaxonomy/$assembly\_geneCalls.fna
    conda deactivate

    conda activate bioinformatics
    kaiju -t $kaiju/nodes.dmp \
          -f $kaiju/nr/kaiju_db_nr.fmi \
          -i contigsTaxonomy/$assembly\_geneCalls.fna \
          -o contigsTaxonomy/$assembly\_geneCalls_kaiju.out \
          -z 10 \
          -v
    kaiju-addTaxonNames -t $kaiju/nodes.dmp \
                        -n $kaiju/names.dmp \
                        -i contigsTaxonomy/$assembly\_geneCalls_kaiju.out \
                        -o contigsTaxonomy/$assembly\_geneCalls_kaiju.names \
                        -r superkingdom,phylum,order,class,family,genus,species
    conda deactivate

    conda activate anvio6.2
    anvi-import-taxonomy-for-genes -i contigsTaxonomy/$assembly\_geneCalls_kaiju.names \
                                   -c anvioDBs/$assembly.db \
                                   -p kaiju \
                                   --just-do-it
    conda deactivate

  else
    echo "Populating genomes complete for" $assembly
  fi
done










####################################################
####################################################
# Run automatic binning algorithms
####################################################
####################################################

screen -S EG_auto_binning
mkdir ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning
mkdir ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2
cd ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2
mapping=~/Everglades/dataEdited/2019_binning/binning_initial/mapping
scripts=/home/GLBRCORG/bpeterson26/Everglades/code/generalUse
scaffolds=~/Everglades/dataEdited/2019_binning/binning_initial/scaffolds
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""
IFS=$'\n'

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  # First need to generate a depth file
  if [ ! -e depth_to_$assembly.txt ]; then
    echo "Summarizing depths for" $assembly
    summarize_command=$(echo -e "jgi_summarize_bam_contig_depths --outputDepth depth_to_$assembly.txt")
    for mappingFile in $(ls $mapping/*_to_$assembly.bam)
    do
      summarize_command=$(echo $summarize_command $mappingFile)
    done
    echo $summarize_command
    eval $summarize_command
  else
    echo "Depth file generation done for" $assembly
  fi
  # Then run the binning
  if [ ! -d $assembly ]; then
    echo "Binning" $assembly
    mkdir $assembly
    metabat2 -i $scaffolds/$assembly\_filtered_assembly.fna \
              -a depth_to_$assembly.txt \
              -o $assembly/metabat_$assembly \
              -m 2000
  else
    echo "Binning of" $assembly "already done"
  fi
done


# Rename bins
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  cd ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2/$assembly
  ls *fa | while read bin
  do
    newBinName=$(echo $bin | sed "s/$assembly./$assembly\_/")
    echo "Renaming" $bin "to" $newBinName
    mv $bin $newBinName
  done
  $scripts/Fasta_to_Scaffolds2Bin.sh -e fa > ../$assembly\_metabat_S2B.tsv
done






####################################################
####################################################
# Generate read profiles
####################################################
####################################################

chmod +x /home/GLBRCORG/bpeterson26/Everglades/code/executables/anvio_read_profiling.sh
condor_submit /home/GLBRCORG/bpeterson26/Everglades/code/submission/anvio_read_profiling_2019.sub





################################################
################################################
# Estimate number of genomes
################################################
################################################

screen -S EG_anvioDBs
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""


assembly=Sed991Mega19
anvi-display-contigs-stats $assembly.db
# Go to this site: http://localhost:8080




################################################
################################################
# Add taxonomic information
################################################
################################################

screen -S EG_scg
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  anvi-run-scg-taxonomy -c $assembly.db
done






################################################
################################################
# Run CONCOCT binning
################################################
################################################

screen -S EG_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
bin_counts=~/Everglades/dataEdited/2019_binning/binning_initial/estimated_number_of_genomes.csv
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs
IFS=$'\n'

# Bin assemblies
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  clusterNumber=`awk -F ',' -v assembly="$assembly" '$1 == assembly { print $3 }' $bin_counts`
  echo -e "Binning "$assembly": "$clusterNumber
  anvi-cluster-contigs -p $assembly.merged/PROFILE.db \
                        -c $assembly.db \
                        -C CONCOCT \
                        -T 16 \
                        --driver concoct \
                        --clusters $clusterNumber \
                        --length-threshold 2000 \
                        --just-do-it
done


################################################
################################################
# Add MetaBat2 info
################################################
################################################

screen -S EG_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs
PYTHONPATH=""
metabat2=~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  anvi-import-collection $metabat2/$assembly\_metabat_S2B.tsv \
                           -c $assembly.db \
                           -p $assembly.merged/PROFILE.db \
                           -C metabat2 \
                           --contigs-mode
   anvi-script-merge-collections -c $assembly.db \
                                   -i $metabat2/$assembly\_metabat_S2B.tsv \
                                   -o $assembly\_collections.tsv
done



################################################
################################################
# Search for hgcA+ bins
################################################
################################################

# First summarize the bins
screen -S EG_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir anvioDB_processing

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  if [ ! -e anvioDB_processing/summaries/$assembly.summary ]; then
    echo "Summarizing binning from" $assembly
    anvi-summarize -c anvioDBs/$assembly.db \
                    -p anvioDBs/$assembly.merged/PROFILE.db \
                    -C CONCOCT \
                    -o anvioDB_processing/$assembly.summary
  else
    echo $assembly "already summarized"
  fi
done

# Get list of original bin names
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDB_processing
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  ls $assembly.summary/bin_by_bin | sed 's/\///' \
      > $assembly\_original_bin_list.txt
done

# Search bins for hgcA
# First we'll want to find the bins that have hgcA.
# We'll have to go through the summaries for that.
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDB_processing
mkdir hgcA_search
rm -f hgcA_search/hgcA_bin_list.txt

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  echo "Looking for hgcA in" $assembly
  cat $assembly\_original_bin_list.txt | while read bin
  do
    if [ -s $assembly.summary/bin_by_bin/$bin/$bin-hgcaAnvio-hmm-sequences.txt ]; then
      echo $assembly$'\t'$bin >> hgcA_search/original_hgcA_bin_list.txt
    fi
  done
done








################################################
################################################
# Manually bin hgcA+ bins
################################################
################################################

screen -S EG_anvioDBs_binning
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
metabat2=~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2
cd ~/Everglades/dataEdited/2019_binning/binning_initial
# Copy the database folder
#cp -avr anvioDBs anvioDBs_modified

assembly=Sed996Meta19
bin=Bin_9
anvi-refine -p anvioDBs_modified/$assembly.merged/PROFILE.db \
            -c anvioDBs_modified/$assembly.db \
            -C CONCOCT \
            -b $bin \
            -A anvioDBs_modified/$assembly\_collections.tsv \
            --taxonomic-level "t_phylum"



################################################
################################################
# Summarize/export curated bins
################################################
################################################

##########################
# Rename and summarize bins
##########################
screen -S EG_2019_binsPostProcessing
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
PERL5LIB=""
IFS=$'\n'
cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs_modified

# Rename them
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  # If you need to delete the renamed collection to make a new one:
  # anvi-delete-collection -p $assembly.db -C refined_and_renamed
  # rm renamed_bins_key.txt
  anvi-rename-bins -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    --collection-to-read CONCOCT \
                    --collection-to-write refined_and_renamed \
                    --prefix $assembly \
                    --report-file renamed_bins_key.txt
done

# Summarize them
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  # If old summary exists that we want to delete, uncomment
  # the following:
  # rm -r $assembly.summary.curated
  if [ ! -d $assembly.curated.summary ]; then
    echo "Summarizing bins for" $assembly
    anvi-summarize -c $assembly.db \
                    -p $assembly.merged/PROFILE.db \
                    -C refined_and_renamed \
                    -o $assembly.curated.summary
  else
    echo "We already summarized the curated bins for" $assembly
  fi
done
conda deactivate


##########################
# Pull out DNA files from hgcA+ bins
##########################
# First need to set up new directory
cd ~/Everglades/dataEdited/2019_binning/binning_initial/
mkdir binsRaw
mkdir binsRaw/DNA
binsRaw=~/Everglades/dataEdited/2019_binning/binning_initial/binsRaw

for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  if [ -e ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs_modified/$assembly.curated.summary ]; then
    if [ ! -e $binsRaw/DNA/$assembly* ]; then
      cd ~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs_modified/$assembly.curated.summary/bin_by_bin
      for bin in $(ls | sed 's/\///')
      do
        isThereHgcA=`cat $bin/$bin\-hgcaAnvio-hmm-sequences.txt | wc -l`
        if [ ! $isThereHgcA -eq 0 ]; then
          echo "Copying" $bin "to binsRaw folder"
          cp $bin/$bin-contigs.fa $binsRaw/DNA/$bin.fna
        else
          echo "No hgcA in" $bin
        fi
      done
    else
      echo "Hey, there are some bins from" $assembly "already in here"
      echo "You might wanna check that out before you start overwriting stuff"
    fi
  else
    echo "Summarize anvioDB for" $assembly", dummy."
  fi
done
