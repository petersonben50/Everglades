
#########################
# Prep fasta files
#########################
# Got hgcA nucleic acid sequences from Elizabeth (from her paper), and uploaded them:
# /home/GLBRCORG/bpeterson26/Everglades/references/hgcA/mcdaniel_2020/hgcA_coding_regions.fna

#########################
# Get set up
#########################
screen -S EG_sortmehgcA
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate seqtk

mkdir /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/hgcA_analysis
mkdir /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/hgcA_analysis/workingDirectory
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/hgcA_analysis/workingDirectory

referenceDB=/home/GLBRCORG/bpeterson26/references/sortmerna_dbs
metagenomeList=/home/GLBRCORG/bpeterson26/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt
linesToCut=1600000
metagenomeFolder=/home/GLBRCORG/bpeterson26/Everglades/dataEdited/metagenomes



#########################
# Interleave and split files
#########################
cat $metagenomeList | while read metagenome
do
  echo "Working on" $metagenome
  seqtk mergepe $metagenomeFolder/$metagenome\_R1.fastq.gz \
                $metagenomeFolder/$metagenome\_R2.fastq.gz | \
  split -a 4 -l $linesToCut - $metagenome\_splitFiles_
done
ls *_splitFiles_* > split_file_list.txt



#########################
# Run sortmerna on metagenomes to identify hgcA
#########################
#chmod +x /home/GLBRCORG/bpeterson26/Everglades/code/executables/sortmerna_featureOfInterest.sh
condor_submit /home/GLBRCORG/bpeterson26/Everglades/code/submission/sortmerna_featureOfInterest_hgcA.sub
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/hgcA_analysis/workingDirectory
mkdir mapping_to_hgcA non_mapping_to_hgcA
mv hgcA_mappingTo.* mapping_to_hgcA
mv hgcA_nonMapping.* non_mapping_to_hgcA
rm -f *_splitFiles_*


#########################
# Process the output
#########################
screen -S EG_hgcA_reads
metagenomeList=/home/GLBRCORG/bpeterson26/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt
cd /home/GLBRCORG/bpeterson26/Everglades/dataEdited/2019_readBasedAnalysis/hgcA_analysis
mkdir sortmerna_output
echo -e "metagenomeID\thgcA_reads\tnon_hgcA_reads" > hgcA_reads_count.txt
cat $metagenomeList | while read metagenome
do
  echo "Working on" $metagenome

  #cat workingDirectory/mapping_to_hgcA/hgcA_mappingTo.$metagenome* > sortmerna_output/$metagenome\_hgcA_reads.fastq
  hgcA_reads=`wc -l sortmerna_output/$metagenome\_hgcA_reads.fastq`

  #cat workingDirectory/non_mapping_to_hgcA/hgcA_nonMapping.$metagenome* > sortmerna_output/$metagenome\_non_hgcA_reads.fastq
  nonHgcA_reads=`wc -l sortmerna_output/$metagenome\_non_hgcA_reads.fastq`

  echo -e $metagenome"\t"$hgcA_reads"\t"$nonHgcA_reads >> hgcA_reads_count.txt
done
