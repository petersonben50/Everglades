#!/bin/sh

######################
# code/assembly_processing.sh
# Benjamin D. Peterson

# This set of scripts will contain everything
# needed for the initial processing of the
# metagenomes.
######################


############################################
############################################
# Data retrieval and inspection
############################################
############################################

######################
# KMBP004: First round of metagenomes
######################

# Get set up
mkdir ~/Everglades/dataRaw/metagenomes


# Now retrieve the metagenomes using lftp
# (in the lftp virtual environment)
# Log onto GLBRC server
#screen -S Everglades_MG_transfer
#cd Everglades/dataRaw/metagenomes

#/opt/wget-1.19/bin/wget -r -np -nH --ask-password ftps://mcmahonk@titan.bch.msu.edu/20180921_DNASeq_PE150
#mv 20180921_DNASeq_PE150/* ./
#rm -fr 20180921_DNASeq_PE150
#exit


# Let's get a list of our metagenomes.
# Save out a list of all the metagenome IDs here:



############################################
############################################
# Trimming metagenomes
############################################
############################################

screen -S Everglades_metagenome_trimming
mkdir ~/Everglades/dataEdited
mkdir ~/Everglades/dataEdited/metagenomes
mkdir ~/Everglades/dataEdited/metagenomes/reports

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""



######################
# Trim metagenomes
######################

cd ~/Everglades/dataRaw/metagenomes

read_storage=~/Everglades/dataEdited/metagenomes
ancillary_info=~/Everglades/dataEdited/metagenomes/reports

cat ~/Everglades/metadata/lists/metagenome_list.csv | while read metagenome
do

  if [ ! -e $read_storage/$metagenome\_R1.fastq.gz ]; then

    echo "Processing" $metagenome

    fastp --in1 $metagenome*R1*fastq.gz \
          --in2 $metagenome*R2*fastq.gz \
          --out1 $read_storage/$metagenome\_R1.fastq.gz \
          --out2 $read_storage/$metagenome\_R2.fastq.gz \
          --unpaired1 $read_storage/$metagenome\_single.fastq.gz \
          --unpaired2 $read_storage/$metagenome\_single.fastq.gz \
          --merge \
          --merged_out $read_storage/$metagenome\_merged.fastq.gz \
          --failed_out $ancillary_info/$metagenome\_failed.fastq.gz \
          --html $ancillary_info/$metagenome\_report.html \
          --detect_adapter_for_pe \
          --cut_tail \
          --cut_tail_window_size 10 \
          --cut_tail_mean_quality 20 \
          --length_required 100


    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_R1.fastq.gz
    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_R2.fastq.gz
    fastqc -o $ancillary_info \
            $read_storage/$metagenome\_single.fastq.gz

  else

    echo "Already processed" $metagenome

  fi

done



















######################
# Count reads in metagenome
######################

screen -S Everglades_metagenome_read_counting

cd ~/Everglades/dataRaw/metagenomes

read_storage=~/Everglades/dataEdited/metagenomes
ancillary_info=~/Everglades/dataEdited/metagenomes/reports

echo -e "metagenomeID\tforwardReads\treverseReads\tsingleReads\tmergedReads" > $ancillary_info/metagenome_read_count.tsv

cat ~/Everglades/metadata/lists/metagenome_list.csv | while read metagenome
do

  echo "Working on" $metagenome

  forwardCount=$(zgrep -c "^@" $read_storage/$metagenome\_R1.fastq.gz)
  reverseCount=$(zgrep -c "^@" $read_storage/$metagenome\_R2.fastq.gz)
  singleCount=$(zgrep -c "^@" $read_storage/$metagenome\_single.fastq.gz)
  mergeCount=$(zgrep -c "^@" $read_storage/$metagenome\_merged.fastq.gz)

  echo -e $metagenome"\t"$forwardCount"\t"$reverseCount"\t"$singleCount"\t"$mergeCount >> $ancillary_info/metagenome_read_count.tsv

done








############################################
############################################
# Metagenome assembly
############################################
############################################

mkdir ~/Everglades/code/assembly_processing

######################
# Assemblies by metaSPades
######################

screen -S Everglades_metagenome_coassembly

cd ~/Everglades/dataEdited/assemblies
mkdir ~/Everglades/dataEdited/assemblies
mkdir ~/Everglades/dataEdited/assemblies/assembly_files/

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

code=~/Everglades/code/assembly_processing/
assembly_grouping=~/Everglades/metadata/lists/assembly_groups_metaspades.csv
read_storage=~/Everglades/dataEdited/metagenomes/
output=~/Everglades/dataEdited/assemblies/assembly_files/


# Make sure the list of wanted assemblies by metaSPADes is up to date
# in ~/Everglades/dataEdited/assemblies/assembly_group_metaspades.csv

# Get a list of group names
cd ~/Everglades/metadata/lists/
tail -n +2 assembly_groups_metaspades.csv | \
    awk -F ',' '{ print $1 }' | \
    uniq \
    > assembly_list_metaspades.txt

cd ~/Everglades/dataEdited/assemblies

cat ~/Everglades/metadata/lists/assembly_list_metaspades.txt | while read assembly
do

  if [ ! -d $output/$assembly ]; then
    mkdir $output/$assembly
  fi

  # Run group assemblies with metaSPADes

  if [ ! -e $output/$assembly/scaffolds.fasta ]; then
    echo "Assembling" $assembly
    python $code/assembly_by_group_metaspades.py $assembly \
                                                  $assembly_grouping \
                                                  $read_storage \
                                                  $output/$assembly
    # To continue a paused run:
    # cd $output
    # metaspades.py --continue -o $assembly

  else
    echo $assembly "already assembled"
  fi

done


#exit
