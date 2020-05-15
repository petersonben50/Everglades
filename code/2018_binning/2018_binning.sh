#!/bin/sh

##########################
# code/2018_binning/2018_binning.sh
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
mkdir ~/Everglades/dataEdited/2018_binning
cd ~/Everglades/dataEdited
mkdir 2018_binning/binning_initial
mkdir 2018_binning/binning_initial/scaffolds
export PATH=/home/GLBRCORG/bpeterson26/miniconda3/bin:$PATH
source activate anvio5
PYTHONPATH=/home/GLBRCORG/bpeterson26/miniconda3/envs/anvio5/lib/python3.6/site-packages/
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt)
do
  if [ ! -e 2018_binning/binning_initial/scaffolds/$assembly\_filtered_assembly.fna ]; then
    echo "Processing" $assembly "scaffolds for binning"
      anvi-script-reformat-fasta assemblies/scaffolds/$assembly\_assembly.fna \
                              -o 2018_binning/binning_initial/scaffolds/$assembly\_filtered_assembly.fna \
                              -l 2000
  else
    echo $assembly "scaffolds already processed for binning"
  fi
done

#exit



##########################
# Map reads to filtered scaffolds
##########################

screen -S EG_binning

cd ~/Everglades/dataEdited
mkdir 2018_binning/binning_initial/mapping
mkdir 2018_binning/binning_initial/mapping/indices

mapping=~/Everglades/dataEdited/2018_binning/binning_initial/mapping

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt)
do

  if [ ! -e $mapping/indices/$assembly.1.bt2 ]; then
    echo "Make index for assembly" $assembly
    /opt/bifxapps/bowtie2-2.2.2/bowtie2-build 2018_binning/binning_initial/scaffolds/$assembly\_filtered_assembly.fna \
                                              $mapping/indices/$assembly
  else
    echo "Already made index for" $assembly
  fi

  for metagenome in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt)
  do
    if [ ! -e $mapping/$metagenome\_to_$assembly.bam ]; then
      echo "Mapping MG" $metagenome "to" $assembly
      /opt/bifxapps/bowtie2-2.2.2/bowtie2 -x $mapping/indices/$assembly \
                                          -1 metagenomes/$metagenome\_R1.fastq.gz \
                                          -2 metagenomes/$metagenome\_R2.fastq.gz \
                                          -U metagenomes/$metagenome\_single.fastq.gz,metagenomes/$metagenome\_merged.fastq.gz \
                                          -p 12 \
                                          -S $mapping/$metagenome\_to_$assembly.sam
      echo "Fixing up" $metagenome
      samtools view $mapping/$metagenome\_to_$assembly.sam -o $mapping/raw.$metagenome\_to_$assembly.bam
      samtools sort -m 10G -@ 12 $mapping/raw.$metagenome\_to_$assembly.bam -o $mapping/$metagenome\_to_$assembly.bam
      samtools index $mapping/$metagenome\_to_$assembly.bam

      rm -f $mapping/*.sam $mapping/raw.*.bam

    else
      echo "Already mapped" $metagenome "to" $assembly
    fi
  done

done









####################################################
####################################################
# Run automatic binning algorithms
####################################################
####################################################

screen -S EG_auto_binning
mkdir ~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning
mkdir ~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning/metabat2
cd ~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning/metabat2
mapping=~/Everglades/dataEdited/2018_binning/binning_initial/mapping
scripts=/home/GLBRCORG/bpeterson26/Everglades/code/generalUse
scaffolds=~/Everglades/dataEdited/2018_binning/binning_initial/scaffolds
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do

  # First need to generate a depth file
  if [ ! -e depth_to_$assembly.txt ]; then
    echo "Summarizing depths for" $assembly
    jgi_summarize_bam_contig_depths --outputDepth depth_to_$assembly.txt \
        $mapping/ENP18_001_002_003_to_$assembly.bam \
        $mapping/ENP18_024_025_to_$assembly.bam \
        $mapping/ENP18_030_032_to_$assembly.bam \
        $mapping/ENP18_048_049_50_to_$assembly.bam \
        $mapping/ENP18_061_to_$assembly.bam
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
cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do

  cd ~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning/metabat2/$assembly

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
cd ~/Everglades/dataEdited/2018_binning/binning_initial
mkdir anvioDBs

for assembly in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt)
do
  if [ ! -f anvioDBs/$assembly.db ]; then
    anvi-gen-contigs-database -f scaffolds/$assembly\_filtered_assembly.fna \
                              -o anvioDBs/$assembly.db \
                              -n $assembly
  else
    echo "Contig database for" $assembly "already exists."
  fi
done

source deactivate





###############################
# Populate contigs DB with HMMs
###############################

screen -S EG_anvioDBs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""
hgcAHMM=~/Everglades/references/hgcaAnvio
cd ~/Everglades/dataEdited/2018_binning/binning_initial
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt)
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
cd ~/Everglades/dataEdited/2018_binning/binning_initial
mkdir contigsTaxonomy/
IFS=$'\n'

for assembly in $(cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt)
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
# Generate read profiles
####################################################
####################################################

screen -S EG_anvioDBs
cd ~/Everglades/dataEdited/2018_binning/binning_initial
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do
  if [ ! -d anvioDBs/$assembly.merged ]; then
    cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
    do
      if [ ! -d anvioDBs/$metagenome\_to_$assembly.profile ]; then
        echo "Profiling reads mapped from:" $metagenome "to" $assembly
        anvi-profile -c anvioDBs/$assembly.db \
                      -i mapping/$metagenome\_to_$assembly.bam \
                      --write-buffer-size 1000 \
                      --min-contig-length 2000 \
                      --num-threads 10 \
                      -o anvioDBs/$metagenome\_to_$assembly.profile
      else
        echo "STOP: Profiling reads mapped to:" $assembly "from" $metagenome "is already done"
      fi
    done
    anvi-merge anvioDBs/*to_$assembly.profile/PROFILE.db \
                -o anvioDBs/$assembly.merged \
                -c anvioDBs/$assembly.db \
                -S $assembly\_merged \
                --skip-hierarchical-clustering
  else
    echo "Merging profiles complete for" $groupName
  fi
done







################################################
################################################
# Estimate number of genomes
################################################
################################################

screen -S EG_anvioDBs
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""


assembly=Pw05Meta18
anvi-display-contigs-stats $assembly.db
# Go to this site: http://localhost:8080




################################################
################################################
# Add taxonomic information
################################################
################################################

screen -S EG_scg
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDBs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=""

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
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
bin_counts=~/Everglades/dataEdited/2018_binning/binning_initial/estimated_number_of_genomes.csv
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDBs

# Bin assemblies
cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do
  clusters=$(awk -F ',' -v assembly="$assembly" '$1 == assembly { print $3 }' $bin_counts)
  echo "Binning" $assembly "into" $clusters "clusters"
  anvi-cluster-contigs -p $assembly.merged/PROFILE.db \
                        -c $assembly.db \
                        -C CONCOCT \
                        -T 16 \
                        --driver concoct \
                        --clusters $clusters \
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
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDBs
PYTHONPATH=""
metabat2=~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning/metabat2

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
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
cd ~/Everglades/dataEdited/2018_binning/binning_initial
mkdir anvioDB_processing

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
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
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDB_processing
cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do
  ls $assembly.summary/bin_by_bin | sed 's/\///' \
      > $assembly\_original_bin_list.txt
done

# Search bins for hgcA
# First we'll want to find the bins that have hgcA.
# We'll have to go through the summaries for that.
cd ~/Everglades/dataEdited/2018_binning/binning_initial/anvioDB_processing
mkdir hgcA_search
rm -f hgcA_search/hgcA_bin_list.txt

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
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
metabat2=~/Everglades/dataEdited/2018_binning/binning_initial/autoBinning/metabat2
cd ~/Everglades/dataEdited/2018_binning/binning_initial
# Copy the
#cp -avr anvioDBs anvioDBs_modified

assembly=Pw03Meta18
bin=Bin_1
anvi-refine -p anvioDBs_modified/$assembly.merged/PROFILE.db \
            -c anvioDBs_modified/$assembly.db \
            -C CONCOCT \
            -b $bin \
            -A anvioDBs_modified/$assembly\_collections.tsv \
            --taxonomic-level "t_phylum"
