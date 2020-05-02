#!/bin/sh

######################
# code/2019_analysis_assembly/assembly_mapping_2019.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Mapping reads to scaffolds for 2019
############################################
############################################

######################
# Map reads and process output
######################

screen -S EG_assembly_mapping

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

mkdir ~/Everglades/dataEdited/2018_analysis_assembly
cd ~/Everglades/dataEdited/2018_analysis_assembly
mkdir mapping
mkdir mapping/indices

read_storage=~/Everglades/dataEdited/metagenomes

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do

  if [ -e ~/Everglades/dataEdited/assemblies/scaffolds/$assembly\_assembly.fna ]; then

    echo $assembly "has been cleaned, let's map to it"

    if [ ! -e mapping/indices/bowtie_index_$assembly.1.bt2 ]; then
      echo "Building index for" $assembly
      bowtie2-build ~/Everglades/dataEdited/assemblies/scaffolds/$assembly\_assembly.fna \
                    mapping/indices/bowtie_index_$assembly
    else
      echo "Already built index for" $assembly
    fi

    cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
    do

      if [ ! -e mapping/$metagenome\_to_$assembly.bam ]; then

        echo "Mapping" $metagenome "to" $assembly
        bowtie2 -x mapping/indices/bowtie_index_$assembly \
                -1 $read_storage/$metagenome\_R1.fastq.gz \
                -2 $read_storage/$metagenome\_R2.fastq.gz \
                -U $read_storage/$metagenome\_merged.fastq.gz,$read_storage/$metagenome\_single.fastq.gz \
                -p 12 \
                -S mapping/$metagenome\_to_$assembly.sam

        samtools view mapping/$metagenome\_to_$assembly.sam \
                      -o mapping/$metagenome\_to_$assembly.unsorted.bam
        # Sort the BAM files
        samtools sort -m 10G \
                      -@ 8 \
                      mapping/$metagenome\_to_$assembly.unsorted.bam \
                      -o mapping/$metagenome\_to_$assembly.bam
        # Index the BAM files.
        samtools index mapping/$metagenome\_to_$assembly.bam

        # Calculate depth at each amino acid residue
        #samtools depth -aa mapping/$metagenome\_to_$assembly.bam > mapping/$metagenome\_to_$assembly.depth

      else
        echo "Mapping of" $metagenome "to" $assembly "is already done"
      fi
    done
  else
    echo "Still gotta assemble and clean" $assembly
  fi
done

cd mapping
rm -f *.sam *unsorted.bam
