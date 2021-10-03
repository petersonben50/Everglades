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
# Run MetaBat2
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


# Rename bins and copy them
mkdir ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2/all_bins

cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  cd ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2/$assembly
  ls *fa | while read bin
  do
    newBinName=$(echo $bin | sed "s/$assembly./$assembly\_/")
    echo "Renaming" $bin "to" $newBinName
    #mv $bin $newBinName
    cp $newBinName ../all_bins/$newBinName
  done
  #$scripts/Fasta_to_Scaffolds2Bin.sh -e fa > ../$assembly\_metabat_S2B.tsv
done


##########################
# Running Prodigal on autogenerated bins
##########################
screen -S EG_bin_ORF_prediction
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

cd ~/Everglades/dataEdited/2019_binning/binning_initial/autoBinning/metabat2/all_bins
mkdir ../ORFs
ls *fa | sed 's/.fa//' | while read bin
do
  echo "Working on" $bin
  if [ ! -e ../ORFs/$bin.gff ]; then
    prodigal -i $bin.fa \
              -o ../ORFs/$bin.gff \
              -f gff \
              -a ../ORFs/$bin.faa \
              -d ../ORFs/$bin.fna \
              -p single
  else
    echo $bin "already run."
  fi
done


cd ../ORFs
cat ../good_bins_list.txt | while read bin
do
  echo "Cleaning" $bin
  python $scripts/cleanFASTA.py $bin.fna
  mv -f $bin.fna_temp.fasta $bin.fna
  python $scripts/cleanFASTA.py $bin.faa
  mv -f $bin.faa_temp.fasta $bin.faa
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

# Generate list of hgcA+ bins
cd $binsRaw/DNA
ls *.fna | \
  sed 's/.fna//' \
  > binsRaw_hgcA_list.txt



####################################################
####################################################
# Check quality of bins
####################################################
####################################################

##########################
# Completeness/redundancy estimates from anvio
##########################
binsRaw=~/Everglades/dataEdited/2019_binning/binning_initial/binsRaw
mkdir $binsRaw/anvio_data
mkdir $binsRaw/anvio_data/completeness_redundancy

# Copy summary files into a single folder.
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  summary=~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs_modified/$assembly.curated.summary
  if [ -e $summary/bins_summary.txt ]; then
    cp $summary/bins_summary.txt $binsRaw/anvio_data/completeness_redundancy/$assembly\_bins_summary.txt
  else
    echo $assembly "has not been summarized."
  fi
done

# Concatenate summaries into a single file.
cd $binsRaw/anvio_data/completeness_redundancy
head -n 1 Sed994Meta19_bins_summary.txt > hgcA_bins_summary_all.txt
for file in $(ls *_bins_summary.txt)
do
  tail -n +2 $file >> bins_summary_all.txt
done

# Only keep the summaries for the hgcA+ bins.
head -n 1 bins_summary_all.txt > bins_summary_hgcA.txt
for hgcA_bin in $(cat $binsRaw/DNA/binsRaw_hgcA_list.txt)
do
  grep $hgcA_bin bins_summary_all.txt >> bins_summary_hgcA.txt
done

head -n 1 hgcA_bins_summary.txt > bins_summary_hgcA_good.txt
awk '{ if (($7 > 50) && ($8 < 10)) print $0 }' bins_summary_hgcA.txt >> bins_summary_hgcA_good.txt
tail -n +2 bins_summary_hgcA_good.txt | \
  awk '{ print $1 }' \
  > bins_list_hgcA_good.txt

# Download bins_summary_hgcA.txt to local computer
# dataEdited/2019_binning/binning_initial/binsRaw


##########################
# Completeness/redundancy estimates from CheckM
##########################

screen -S EG_checkM
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics

cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsRaw
if [ -d checkM ]; then
  echo "Removing old checkM folder"
  rm -rf checkM
fi
mkdir checkM
checkm lineage_wf \
        -x .fna \
        -t 16 \
        DNA \
        checkM
checkm qa checkM/lineage.ms \
          checkM \
          -o 2 \
          -f checkM/checkm.out \
          --tab_table
awk -F '\t' \
    -v OFS=',' \
    '{ print $1,$6,$7,$8,$9,$11,$13,$15,$17,$19,$23 }' \
    checkM/checkm.out \
    > checkM/checkM_stats.csv

# Download checkM/checkM_stats.csv to local computer:
# dataEdited/2017_analysis_bins/binning/rawBins/bin_quality
cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir binsGood
mkdir binsGood/DNA
mkdir binsGood/checkM
awk -F ',' '{ if (($2 > 50) && ($3 < 10)) print $0 }' \
    binsRaw/checkM/checkM_stats.csv \
    > binsGood/checkM/good_bins_data.txt
awk -F ',' '{ print $1 }' binsGood/checkM/good_bins_data.txt \
  > binsGood/checkM/good_bins_list.txt
cat binsGood/checkM/good_bins_list.txt | while read binsGood
do
  echo "Copying" $binsGood
  cp binsRaw/DNA/$binsGood.fna binsGood/DNA
done

cd binsGood/DNA
scripts=~/Everglades/code/generalUse
ls *fna | while read fna
do
  echo "Cleaning" $fna
  python $scripts/cleanFASTA.py $fna
  mv -f $fna\_temp.fasta $fna
done



####################################################
####################################################
# Check out taxonomy of bins with GTDB
####################################################
####################################################

screen -S EG_GTDB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate gtdbtk
cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood
rm -rf taxonomy
mkdir taxonomy

gtdbtk classify_wf \
        --cpus 16 \
        --extension fna \
        --genome_dir ./DNA \
        --out_dir taxonomy
# Summarize them
cd taxonomy
grep -h '_Bin_' gtdbtk.*.summary.tsv \
        | awk -F '\t' '{ print $1"\t"$2 }' \
        > taxonomy_summary.txt
conda deactivate



####################################################
####################################################
# Compare bins by differential coverage
####################################################
####################################################

# Coverage data from anvio
binsGood=~/Everglades/dataEdited/2019_binning/binning_initial/binsGood
mkdir $binsGood/coverage
IFS=$'\n'
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  summary=~/Everglades/dataEdited/2019_binning/binning_initial/anvioDBs_modified/$assembly.curated.summary/bins_across_samples
  if [ -e $summary/mean_coverage_Q2Q3.txt ]; then
    cp $summary/mean_coverage_Q2Q3.txt $binsGood/coverage/$assembly\_coverage.txt
  else
    echo $assembly "has not been summarized."
  fi
done
cd $binsGood/coverage/
head -n 1 Sed991Mega19_coverage.txt > coverage.txt
for assembly in $(cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt)
do
  tail -n +2 $assembly\_coverage.txt >> coverage.txt
done
cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood
head -n 1 coverage/coverage.txt > coverage/coverage_goodBins.txt
grep -f checkM/good_bins_list.txt coverage/coverage.txt >> coverage/coverage_goodBins.txt
# Download the "coverage_goodBins.txt" to my computer




####################################################
####################################################
# Get ORFs for bins
####################################################
####################################################

screen -S EG_binsORFS
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
scripts=/home/GLBRCORG/bpeterson26/Everglades/code/generalUse

cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood
mkdir ORFs
cat good_bins_list.txt | while read bin
do
  if [ ! -e ORFs/$bin.gff ]; then
    prodigal -i DNA/$bin.fna \
              -o ORFs/$bin.gff \
              -f gff \
              -a ORFs/$bin.faa \
              -d ORFs/$bin.fna \
              -p single
  else
    echo $bin "already run."
  fi
done
cd ORFs
cat ../good_bins_list.txt | while read bin
do
  echo "Cleaning" $bin
  python $scripts/cleanFASTA.py $bin.fna
  mv -f $bin.fna_temp.fasta $bin.fna
  python $scripts/cleanFASTA.py $bin.faa
  mv -f $bin.faa_temp.fasta $bin.faa
done



####################################################
####################################################
# Run ANI comparisons on good bins
####################################################
####################################################
# The following workflow is the code to run
# Sarah Stevens's ANI calculator on a
# folder full of bins.
# Details: https://github.com/sstevens2/ani_compare_dag


# Switch to CHTC
mkdir ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison
cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison
wget https://ani.jgi-psf.org/download_files/ANIcalculator_v1.tgz
tar -xzvf ANIcalculator_v1.tgz
git clone https://github.com/sstevens2/ani_compare_dag.git
mv ani_compare_dag EG_bins_ANI
cd EG_bins_ANI/
mkdir goodBins
cp ~/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ORFs/*fna goodBins/
# from GLBRC to CHTC into ANI_comparison/ani_compare_dag/
echo 'goodBins' > groupslist.txt
# Change path of executable and
# transfer_input_files lines
#executable = /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison/EG_bins_ANI/group.sh
#transfer_input_files = /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison/ANIcalculator_v1/ANIcalculator,/home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_binning/binning_initial/binsGood/ANI_comparison/ANIcalculator_v1/nsimscan,$(spllist),$(totransfer)
condor_submit_dag runAllANIcompare.dag
# Download output file (goodBins.all.ani.out.cleaned)
# to my computer:
# dataEdited/2019_binning/binning_initial/binsGood/goodBins.all.ani.out.cleaned

cd ~/Everglades/dataEdited/2019_binning/binning_initial
mkdir binsFinal
mkdir binsFinal/DNA
mkdir binsFinal/ORFs
cat binsFinal_list.txt | while read bin
do
  cp binsGood/DNA/$bin.fna binsFinal/DNA
  cp binsGood/ORFs/$bin* binsFinal/ORFs
done

cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal
scripts=~/Everglades/code/generalUse/
$scripts/Fasta_to_Scaffolds2Bin.sh -e fna \
                                    -i DNA \
                                    > binsFinal_S2B.tsv
cat DNA/*.fna > DNA.fna
$scripts/Fasta_to_Scaffolds2Bin.sh -e faa \
                                    -i ORFs \
                                    > binsFinal_G2B.tsv
cat ORFs/*.faa > ORFs.faa



####################################################
####################################################
# Confirm hgcA in bins
####################################################
####################################################


##########################
# hgcA search
##########################
cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/Everglades/code/generalUse/

hmmsearch --tblout hgcA/hgcA.out \
          --cpu 4 \
          --cut_tc \
          ~/references/hgcA/hgcA.hmm \
          ORFs.faa \
          > hgcA/hgcA_report.txt
grep 'Sed991Mega19_000001090714_3' binsFinal_G2B.tsv
python $scripts/extract_protein_hitting_HMM.py \
        hgcA/hgcA.out \
        ORFs.faa \
        hgcA/hgcA.faa
hmmalign -o hgcA/hgcA.sto \
            ~/references/hgcA/hgcA.hmm \
            hgcA/hgcA.faa
$scripts/convert_stockhold_to_fasta.py hgcA/hgcA.sto

echo -e "hgcA\tbinID" > hgcA/hgcA_bin_key.tsv
grep '>' hgcA/hgcA.faa | \
  sed 's/>//' | \
  while read hgcA
  do
    binID=`awk -F '\t' -v hgcA="$hgcA" '$1 == hgcA {print $2}' binsFinal_G2B.tsv`
    if [ ! -z $binID ]; then
      echo $hgcA "is binned into" $binID
      echo -e "$hgcA\t$binID" >> hgcA/hgcA_bin_key.tsv
    else
      echo $hgcA "is not binned"
    fi
  done



####################################################
####################################################
# Generate phylogenetic tree of bins
####################################################
####################################################




####################################################
####################################################
# Run metabolic HMMs on bins
####################################################
####################################################

screen -S EG_metabolic_HMMs
cd ~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal
scripts=~/Everglades/code/metagenomeBinAnalysis
metabolic_HMMs=~/Everglades/references/metabolic_HMMs
ORFs=~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal/ORFs.faa
workingDirectory=~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=''
PERL5LIB=''
chmod +x $scripts/batch_HMMs.py

python $scripts/batch_HMMs.py --orf_file $workingDirectory/ORFs.faa \
                              --g2b $workingDirectory/binsFinal_G2B.tsv \
                              --hmm_folder $metabolic_HMMs\
                              --hmm_csv $metabolic_HMMs.csv \
                              --output ~/Everglades/dataEdited/2019_binning/batch_HMM_output
conda deactivate
exit


# Look for MHCs
cd $workingDirectory
mkdir MHCs
#chmod +x ~/Everglades/code/metagenomeBinAnalysis/Find_multiheme_protein.py
~/Everglades/code/metagenomeBinAnalysis/Find_multiheme_protein.py ORFs.faa 3
mv ORFs_3_heme* MHCs

echo -e "binID\tgeneID\themeCount" > MHCs/heme_count_bins.tsv
tail -n +2 MHCs/ORFs_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' binsFinal_G2B.tsv`
  geneCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' MHCs/ORFs_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$geneCount
  echo -e $binID"\t"$geneID"\t"$geneCount >> MHCs/heme_count_bins.tsv
done
