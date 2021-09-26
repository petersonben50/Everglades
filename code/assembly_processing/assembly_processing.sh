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


######################
# KMBP006: Next round of porewater metagenomes
######################

# Let's get a list of our metagenomes.
# Save out a list of all the metagenome IDs here:
#screen -S Everglades_MG_transfer
#cd Everglades/dataRaw/metagenomes
#source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
#conda activate lftp
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n200712_Peterson,rohhei3bohZ2Oov -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'

#mv laneBarcode.html KMBP006_laneBarcode.html
#mv md5sum.txt KMBP006_md5sum.txt

######################
# KMBP005: Sediment metagenomes
######################

# Get set up
cd ~/Everglades/dataRaw/metagenomes

# Now retrieve the metagenomes using lftp
# (in the lftp virtual environment)
screen -S EG_metagenome_retrieval
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate lftp

# Download first attempt at sequencing for KMBP004.
#lftp -c 'set ssl:verify-certificate no set ftp:ssl-protect-data true set ftp:ssl-force true; open -u n200805_150PE_NVS1A_S4_L23_Peterson,zoe3kah4Qua2aim -e "mirror -c; quit" ftp://gslanalyzer.qb3.berkeley.edu:990'
#mv laneBarcode.html KMBP005_laneBarcode.html
#mv md5sum.txt KMBP005_md5sum.txt

cd ~/Everglades/dataRaw/metagenomes
mkdir misnamed_files
mv KMBP_005* misnamed_files

cd misnamed_files
ls KMBP_005* | while read misnamed_file
do
  file_name=$(echo $misnamed_file | sed 's/KMBP_00/KMBP00/g')
  echo "Renaming" $misnamed_file "to" $file_name
  cp $misnamed_file ../$file_name
done


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

cat ~/Everglades/metadata/lists/metagenome_list.txt | while read metagenome
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








############################################
############################################
# Check size of metagenomes
############################################
############################################

######################
# Count reads in metagenome pre-trimming
######################

screen -S Everglades_metagenome_read_counting

cd ~/Everglades/dataRaw/metagenomes

read_storage=~/Everglades/dataRaw/metagenomes
ancillary_info=~/Everglades/dataEdited/metagenomes/reports

echo -e "metagenomeID\tforwardReads\treverseReads" > $ancillary_info/metagenome_read_count_pre_trimming.tsv
cat ~/Everglades/metadata/lists/metagenome_list.csv | while read metagenome
do
  echo "Working on" $metagenome
  forwardCount=$(zgrep -c "^@" $read_storage/$metagenome*_R1*.fastq.gz)
  reverseCount=$(zgrep -c "^@" $read_storage/$metagenome*_R2*.fastq.gz)
  echo -e $metagenome"\t"$forwardCount"\t"$reverseCount >> $ancillary_info/metagenome_read_count_pre_trimming.tsv
done


######################
# Count reads in metagenome post-trimming
######################

screen -S Everglades_metagenome_read_counting_post
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


######################
# Coverage pre-trimming
######################

screen -S EG_MG_coverage_counting_pre
cd ~/Everglades/dataEdited/metagenomes/reports
code=~/Everglades/code/generalUse/readfq-master
read_storage=~/Everglades/dataRaw/metagenomes
IFS=$'\n'

echo -e "metagenomeID\tR1\tR2" > metagenome_coverage_pre_trimming.tsv
for metagenome in $(cat ~/Everglades/metadata/lists/metagenome_list.csv)
do
  echo "Counting coverage in" $metagenome
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenome*_R1*.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenome*_R2*.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $metagenome"\t"$R1_count"\t"$R2_count >> metagenome_coverage_pre_trimming.tsv
done



######################
# Coverage post-trimming
######################

screen -S EG_MG_coverage_counting_post
cd ~/Everglades/dataEdited/metagenomes/reports
code=~/Everglades/code/generalUse/readfq-master
read_storage=~/Everglades/dataEdited/metagenomes
IFS=$'\n'

echo -e "metagenomeID\tR1\tR2\tsingle\tmerged" > metagenome_coverage.tsv
for metagenome in $(cat ~/Everglades/metadata/lists/metagenome_list.csv)
do
  echo "Counting coverage in" $metagenome
  R1_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R1.fastq.gz | \
                awk -F " " '{ print $5 }')
  R2_count=$($code/kseq_fastq_base $read_storage/$metagenome\_R2.fastq.gz | \
                awk -F " " '{ print $5 }')
  single_count=$($code/kseq_fastq_base $read_storage/$metagenome\_single.fastq.gz | \
                awk -F " " '{ print $5 }')
  merged_count=$($code/kseq_fastq_base $read_storage/$metagenome\_merged.fastq.gz | \
                awk -F " " '{ print $5 }')
  echo -e $metagenome"\t"$R1_count"\t"$R2_count"\t"$single_count"\t"$merged_count >> metagenome_coverage.tsv
done




######################
# Generate mash sketches for metagenomes
######################

screen -S EG_mash
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics


# Do this for each metagenome individually

mkdir ~/Everglades/dataEdited/mash_data
mkdir ~/Everglades/dataEdited/mash_data/temp_MG_files
mkdir ~/Everglades/dataEdited/mash_data/sketch_files
cd ~/Everglades/dataEdited/mash_data
read_storage=~/Everglades/dataEdited/metagenomes


cat ~/Everglades/metadata/lists/metagenome_list.csv | while read metagenome
do
  if [ ! -e sketch_files/$metagenome.msh ]; then
    cat $read_storage/$metagenome*.fastq.gz > temp_MG_files/$metagenome.fastq.gz
    mash sketch -S 50 \
                -r \
                -m 2 \
                -k 21 \
                -s 100000 \
                -o sketch_files/$metagenome \
                temp_MG_files/$metagenome.fastq.gz
    rm -f temp_MG_files/$metagenome.fastq.gz
  else
    echo "Already sketched" $metagenome
  fi
done


# Combine sketches
cd ~/Everglades/dataEdited/mash_data
mash paste sketch_files/EG_MG_2019_seds sketch_files/KMBP005*.msh

# Run distance calcuation
#mash dist -S 50 \
#          sketch_files/EG_MG_2019.msh \
#          sketch_files/EG_MG_2019.msh \
#          > EG_MG_2019.dist
mash dist -S 50 \
          sketch_files/EG_MG_2019_seds.msh \
          sketch_files/EG_MG_2019_seds.msh \
          > EG_MG_2019_seds.dist




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









######################
# Assemblies by MegaHit
######################

screen -S Everglades_metagenome_megahit

cd ~/Everglades/dataEdited/assemblies

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

code=~/Everglades/code/assembly_processing/
assembly_grouping=~/Everglades/metadata/lists/assembly_groups_megahit.csv
read_storage=~/Everglades/dataEdited/metagenomes/
output=~/Everglades/dataEdited/assemblies/assembly_files/


# Make sure the list of wanted assemblies by MegaHit is up to date
# in ~/Everglades/metadata/lists/assembly_groups_megahit.csv

# Get a list of group names
cd ~/Everglades/metadata/lists
tail -n +2 assembly_groups_megahit.csv | \
    awk -F ',' '{ print $1 }' | \
    uniq \
    > assembly_list_megahit.txt

# Start assembly
cd $output
cat ~/Everglades/metadata/lists/assembly_list_megahit.txt | while read assembly
do

  if [ ! -e $output/$assembly/final.contigs.fa ]; then
    echo "Assembling" $assembly
    python $code/assembly_by_group_megahit.py $assembly \
                                              $assembly_grouping \
                                              $read_storage \
                                              $output/$assembly

  else
    echo $assembly "already assembled"
  fi
done







############################################
############################################
# Clean up metagenome assemblies, get stats
############################################
############################################

screen -S EG_clean_metagenome_assembly

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate anvio6.2
PYTHONPATH=/home/GLBRCORG/bpeterson26/miniconda3/envs/anvio6.2/lib/python3.6/site-packages/

cd ~/Everglades/dataEdited/assemblies
mkdir scaffolds
mkdir scaffolds/renaming_reports


######################
# Clean Megahit asssemblies
######################

cat ~/Everglades/metadata/lists/assembly_list_megahit.txt | while read assembly
do

  if [ ! -e assembly_files/$assembly/final.contigs.fa ]; then
    echo "Gotta run assembly" $assembly "first, dummy"
  else
    if [ -e scaffolds/$assembly\_assembly.fna ]; then
      echo "Fuck dude, you already cleaned the assembly for" $assembly". Relax."
    else
      echo "Cleaning up assembly" $assembly
      anvi-script-reformat-fasta assembly_files/$assembly/final.contigs.fa \
                                  -o scaffolds/$assembly\_assembly.fna \
                                  -l 1000 \
                                  --simplify-names \
                                  --prefix $assembly \
                                  --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
    fi
  fi
done


######################
# Clean metaSPADes asssemblies
######################

cat ~/Everglades/metadata/lists/assembly_list_metaspades.txt | while read assembly
do

  if [ ! -e assembly_files/$assembly/scaffolds.fasta ]; then
    echo "Gotta run assembly" $assembly "first, dummy"
  else
    if [ -e scaffolds/$assembly\_assembly.fna ]; then
      echo "Fuck dude, you already cleaned the assembly for" $assembly". Relax."
    else
      echo "Cleaning up assembly" $assembly
      anvi-script-reformat-fasta assembly_files/$assembly/scaffolds.fasta \
                                  -o scaffolds/$assembly\_assembly.fna \
                                  -l 1000 \
                                  --simplify-names \
                                  --prefix $assembly \
                                  --report-file scaffolds/renaming_reports/$assembly\_report_file.txt
    fi
  fi
done





######################
# Combine lists of assemblies
######################
cd ~/Everglades/metadata/lists
cat assembly_list_megahit.txt \
    assembly_list_metaspades.txt \
    > assembly_list_all.txt





######################
# Get assembly stats
######################

screen -S EG_assembly_stats

cd ~/Everglades/dataEdited/assemblies
mkdir reports
scripts=~/Everglades/code/assembly_processing

# Get stats on all assemblies
cat ~/Everglades/metadata/lists/assembly_list_all.txt | while read assembly
do
  if [ -e scaffolds/$assembly\_assembly.fna ]; then
    if [ ! -e reports/$assembly\_report.txt ]; then
      echo $assembly "has been cleaned. Let's get some stats on it."
      perl $scripts/abyss-fac.pl scaffolds/$assembly\_assembly.fna \
          > reports/$assembly\_report.txt
    else
      echo "Already checked" $assembly
    fi
  else
    echo $assembly "has not been cleaned. Do that first."
  fi
done

# Aggregate stats
cd ~/Everglades/dataEdited/assemblies/reports
echo -e "n\tn:200\tL50\tmin\tN80\tN50\tN20\tmax\tsum\tassemblyID" > all_assemblies_stats.txt
for file in *report.txt
do
  tail -n +2 $file >> all_assemblies_stats.txt
done

# Clean up the report file
cat ~/Everglades/metadata/lists/assembly_list_all.txt | while read assembly
do
  sed "s/scaffolds\/$assembly\_assembly.fna/$assembly/" all_assemblies_stats.txt \
    > all_assemblies_stats.txt_edited
  mv -f all_assemblies_stats.txt_edited all_assemblies_stats.txt
done

# exit







############################################
############################################
# Predict open reading frames
############################################
############################################

screen -S EG_predict_ORFs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

cd ~/Everglades/dataEdited/assemblies
mkdir ORFs

cat ~/Everglades/metadata/lists/assembly_list_all.txt | while read assembly
do
  if [ -e scaffolds/$assembly\_assembly.fna ]; then
    echo $assembly "has been cleaned. Let's predict ORFs."
    if [ -e ORFs/$assembly.faa ]; then
      echo "You already predicted ORFs for" $assembly". Relax."
    else
      echo "Time to predict some proteins for" $assembly
      prodigal -i scaffolds/$assembly\_assembly.fna \
                -o ORFs/$assembly.gff \
                -f gff \
                -a ORFs/$assembly.faa \
                -d ORFs/$assembly.fna \
                -p meta

      python ~/Everglades/code/generalUse/cleanFASTA.py ORFs/$assembly.fna
      mv -f ORFs/$assembly.fna_temp.fasta ORFs/$assembly.fna
      python ~/Everglades/code/generalUse/cleanFASTA.py ORFs/$assembly.faa
      mv -f ORFs/$assembly.faa_temp.fasta ORFs/$assembly.faa
    fi
  else
    echo "You gotta clean this shit:" $assembly
  fi
done



######################
# Count ORFs
######################

screen -S EG_count_ORFs

cd ~/Everglades/dataEdited/assemblies

echo -e 'assemblyID\tORF_count' > reports/ORF_counts.tsv
cat ~/Everglades/metadata/lists/assembly_list_all.txt | while read assembly
do
  if [ -e ORFs/$assembly.faa ]; then
    echo "Count ORFs in" $assembly
    ORF_count=$(grep '>' ORFs/$assembly.faa | wc -l)
    echo -e "$assembly\t$ORF_count" >> reports/ORF_counts.tsv
  fi
done
