#!/bin/sh

######################
# code/2019_analysis_assembly/mapping_assembly.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Mapping reads to scaffolds for 2019
############################################
############################################
cd ~/Everglades/metadata/lists/
grep '19' assembly_list_all.txt > 2019_analysis_assembly_list.txt
grep 'KMBP00[5-6]' metagenome_list.csv > 2019_analysis_assembly_metagenomes_list.txt

######################
# Map reads and process output
######################

screen -S EG_assembly_mapping

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

mkdir ~/Everglades/dataEdited/2019_analysis_assembly
cd ~/Everglades/dataEdited/2019_analysis_assembly
mkdir mapping
mkdir mapping/indices

read_storage=~/Everglades/dataEdited/metagenomes

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
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

    cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
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



######################
# Calculate fraction of reads mapped to each assembly
######################
mkdir reports
echo -e "metagenomeID\tassemblyID\tpercentageMapped" > reports/mapping_percentages.tsv
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
  do
    if [ ! -e reports/$metagenome\_to_$assembly.stats ]; then
      echo "Calculating flagstat for mapping of" $metagenome "to" $assembly
      samtools flagstat $metagenome\_to_$assembly.bam > reports/$metagenome\_to_$assembly.stats
    fi
    echo "Aggregating percentrage" $metagenome "to" $assembly
    percentMapped=$(cat reports/$metagenome\_to_$assembly.stats | \
      grep "mapped (" | \
      awk -F "[(|%]" '{print $2}')
    echo -e $metagenome"\t"$assembly"\t"$percentMapped >> reports/mapping_percentages.tsv
  done
done


######################
# Calculate overall fraction of reads in a metagenome mapping to any assembly
######################

screen -S EG_assembly_mapping

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PERL5LIB=""

cd ~/Everglades/dataEdited/2019_analysis_assembly/mapping
mkdir mapped_reads


cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
do
  cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
  do
    if [ ! -e  mapped_reads/$metagenome\_to_$assembly\_singleReads.txt ]; then
      echo "Extracting single reads mapped from" $metagenome "to" $assembly
      samtools view -F 5 $metagenome\_to_$assembly.bam | \
        awk -F '\t' '{ print $1 }' \
        > mapped_reads/$metagenome\_to_$assembly\_singleReads.txt
    fi

    if [ ! -e  mapped_reads/$metagenome\_to_$assembly\_forwardReads.txt ]; then
      echo "Extracting forward reads mapped from" $metagenome "to" $assembly
      samtools view -f 65 -F 4 $metagenome\_to_$assembly.bam | \
        awk -F '\t' '{ print $1 }' \
        > mapped_reads/$metagenome\_to_$assembly\_forwardReads.txt
    fi

    if [ ! -e  mapped_reads/$metagenome\_to_$assembly\_reverseReads.txt ]; then
      echo "Extracting reverse reads mapped from" $metagenome "to" $assembly
      samtools view -f 129 -F 4 $metagenome\_to_$assembly.bam | \
        awk -F '\t' '{ print $1 }' \
        > mapped_reads/$metagenome\_to_$assembly\_reverseReads.txt
    fi

    if [ ! -e  mapped_reads/$metagenome\_to_$assembly\_unmappedReads.txt ]; then
      echo "Extracting reads from" $metagenome "that did not map to" $assembly
      samtools view -f 4 $metagenome\_to_$assembly.bam | \
        awk -F '\t' '{ print $1 }' \
        > mapped_reads/$metagenome\_to_$assembly\_unmappedReads.txt
    fi
  done
done

# Concatenate unique reads
cd ~/Everglades/dataEdited/2019_analysis_assembly/mapping/mapped_reads
cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
do
  echo "Forward reads: "$metagenome
  cat $metagenome\_*forwardReads.txt | \
    sort | uniq \
    > $metagenome\_forwardReads.txt
  echo "Reverse reads: "$metagenome
  cat $metagenome\_*reverseReads.txt | \
    sort | uniq \
    > $metagenome\_reverseReads.txt
  echo "Single reads: "$metagenome
  cat $metagenome\_*singleReads.txt | \
    sort | uniq \
    > $metagenome\_singleReads.txt
done

# Count the number of lines in each type of file
cd ~/Everglades/dataEdited/2019_analysis_assembly/mapping/mapped_reads
rm -rf ../uniq_mapped_reads_MG.tsv
cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
do
  forwardReadCounts=`cat $metagenome\_forwardReads.txt | wc -l`
  reverseReadCounts=`cat $metagenome\_reverseReads.txt| wc -l`
  singleReadCounts=`cat $metagenome\_singleReads.txt | wc -l`
  echo -e $metagenome"\t"$forwardReadCounts"\t"$reverseReadCounts"\t"$singleReadCounts >> ../reports/uniq_mapped_reads_MG.tsv
done



# Clean up
cd ~/Everglades/dataEdited/2019_analysis_assembly/mapping
rm -rf mapped_reads
