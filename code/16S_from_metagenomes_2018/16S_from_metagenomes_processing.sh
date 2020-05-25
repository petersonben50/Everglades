#!/bin/sh

######################
# code/16S_from_metagenomes_processing.sh
# Benjamin D. Peterson

# This set of scripts will assemble 16S
# sequences from all the Hells Canyon
# metagenomes, cluster them, and assign
# taxonomy.
######################

mkdir ~/Everglades/dataEdited/16S_from_MG_2018
mkdir ~/Everglades/dataEdited/16S_from_MG_2018/setup
mkdir ~/Everglades/dataEdited/16S_from_MG_2018/emirge_output



##########################################################
##########################################################
# Assemble needed sequences using EMIRGE
##########################################################
##########################################################


######################
# Run EMIRGE on all samples
######################

screen -S EG_emirge

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate emirge
PYTHONPATH=''

read_storage=~/Everglades/dataEdited/metagenomes
emirge_db=~/references/emirge_db
emirge_output=~/Everglades/dataEdited/16S_from_MG_2018/emirge_output

cd $emirge_output

cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do
  if [ ! -d $emirge_output/$metagenome ]; then
    echo "Running EMIRGE on" $metagenome
    gunzip < $read_storage/$metagenome\_R1.fastq.gz \
        > $emirge_output/$metagenome\_R1.fastq
    gunzip < $read_storage/$metagenome\_R2.fastq.gz \
        > $emirge_output/$metagenome\_R2.fastq
    insert_size=$(awk -F '\t' -v metagenome="$metagenome" '$1 == metagenome { print $2}' ~/Everglades/dataEdited/16S_from_MG_2018/setup/MG_insert_sizes_PW.txt)
    insert_sd=$(awk -F '\t' -v metagenome="$metagenome" '$1 == metagenome { print $3}' ~/Everglades/dataEdited/16S_from_MG_2018/setup/MG_insert_sizes_PW.txt)

    emirge.py $emirge_output/$metagenome \
              -1 $emirge_output/$metagenome\_R1.fastq \
              -2 $emirge_output/$metagenome\_R2.fastq \
              -f ~/references/emirge_db/SILVA_138_SSURef_NR99_tax_silva_trunc.ge1200bp.le2000bp.0.97.fixed.fasta \
              -b ~/references/emirge_db/SILVA_138_SSURef_NR99_tax_silva_trunc.ge1200bp.le2000bp.0.97.fixed \
              -l 150 \
              -i $insert_size \
              -s $insert_sd \
              -a 16 \
              --phred33
    rm -f $emirge_output/$metagenome\_R1.fastq \
          $emirge_output/$metagenome\_R2.fastq
  else
    echo "EMIRGE was already run on" $metagenome
  fi
done











######################
# Pull out all fasta sequences
######################

screen -S EG_emirge_processing

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate emirge
PYTHONPATH=''

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_output
mkdir ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data
mkdir ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/fasta_seqs
output=~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/fasta_seqs

# Gather all fasta sequences
cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do
  cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_output/$metagenome
  echo $metagenome
  # Run rename script to included prior abundances and rank sequences according to abundance
  emirge_rename_fasta.py iter.40 > $metagenome.fasta

  python ~/Everglades/code/generalUse/cleanFASTA.py $metagenome.fasta

  # Add metagenome ID to ensure seq names are unique
  sed "s/>/>$metagenome\_/" $metagenome\_temp.fasta > $output/$metagenome.fasta
done
cd $output
cat *fasta > 16S_seqs_all.fna





######################
# Cluster 16S sequences across metagenomes
######################

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/fasta_seqs

cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i 16S_seqs_all.fna \
              -o 16S_seqs_derep.fna \
              -c 0.97 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl 16S_seqs_derep.fna.clstr > 16S_seqs_cluster.tsv

rm -f 16S_seqs_key.tsv
# Generate a key to rename sequences that cluster together
tail -n +2 16S_seqs_cluster.tsv | awk -F '\t' '{ print $1 }' | while read seqID
do

  clusterID=$(awk -F '\t' -v seqID="$seqID" '$1 == seqID { print $2; exit }' 16S_seqs_cluster.tsv)

  repID=$(awk -F '\t' -v clusterID="$clusterID" '$2 == clusterID && $5 == 1 { print $1 }' 16S_seqs_cluster.tsv)

  echo -e $seqID'\t'$repID
  echo -e $seqID'\t'$repID >> 16S_seqs_key.tsv

done









######################
# Extract abundance of 16S sequences
######################

screen -S EG_emirge_processing
cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/fasta_seqs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
# Extract fraction of reads mapping to each 16S sequence,
# as well as the total number of reads mapping to 16S sequences
# in general.

mkdir ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/coverage
coverage_output=~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/coverage

echo -e 'metagenomeID\ttotal16Sreads' > $coverage_output/mapped_reads.txt
echo -e 'metagenomeID\tseqID\trawAbundance\tlength\tabundance' > $coverage_output/fractional_abundance.txt

cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_output/$metagenome
  echo $metagenome

  # Extract mapped number of single reads
  samtools flagstat iter.40/bowtie.iter.40.PE.bam | \
    sed -n 7p | \
    awk -v pwd="${PWD##*/}" '{print pwd, $1}' \
    >> $coverage_output/mapped_reads.txt

  # Extract the sequences names per sample and store in long format .txt file
  grep ">" $metagenome.fasta | \
    sed "s/>/>$metagenome\_/" | \
    awk -v pwd="${PWD##*/}" '{print pwd, $1, $2, $3, $4}' \
    >> $coverage_output/fractional_abundance.txt
done

# Need to multiply priors by total number of reads to get a count file. Probably could just leave it at abundance though, since we don't know the total number of reads mapping to 16S.

# Remove all 'Prior=', 'NormPrior' and 'Length=' labels from the file
sed -i 's/NormPrior=//g' $coverage_output/fractional_abundance.txt
sed -i 's/Prior=//g' $coverage_output/fractional_abundance.txt
sed -i 's/Length=//g' $coverage_output/fractional_abundance.txt
sed -i 's/>//g' $coverage_output/fractional_abundance.txt
sed -i 's/ /\t/g' $coverage_output/fractional_abundance.txt







######################
# Rename seqs in abundance file with representative ID
######################

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/fasta_seqs

rm -f $coverage_output/fractional_abundance_clustered.txt
cp $coverage_output/fractional_abundance.txt $coverage_output/fractional_abundance_clustered.txt

awk -F '\t' '{ print $1 }' 16S_seqs_key.tsv | while read seqID
do

  clusterID=$(awk -F '\t' -v seqID="$seqID" '$1 == seqID { print $2; exit }' 16S_seqs_key.tsv)
  echo "Replacing" $seqID "with" $clusterID

  sed -i "s/$seqID/$clusterID/" $coverage_output/fractional_abundance_clustered.txt

done










##########################################################
##########################################################
# Classify 16S sequences using mothur
##########################################################
##########################################################

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data
mkdir taxonomy
mkdir taxonomy_processing

# Bring in TaxAss scripts
cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/taxonomy_processing
cp ~/programs/16S/taxAssScripts.zip ./
unzip taxAssScripts.zip
chmod +x *.R *.py
cp ~/programs/16S/programsFor16S.zip ./
unzip programsFor16S.zip
tar -zxvf mothur.tar.gz

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
PYTHONPATH=""
conda activate bioinformatics

# Blast OTUs against FW blast database
references=~/references/16S/
blastn -query fasta_seqs/16S_seqs_derep.fna \
        -task megablast \
        -db $references/FW.fasta.db \
        -out taxonomy_processing/OTU.custom.blast \
        -outfmt 11 \
        -max_target_seqs 5

# Reformat the blast results
blast_formatter -archive taxonomy_processing/OTU.custom.blast \
                -outfmt "6 qseqid pident length qlen qstart qend" \
                -out taxonomy_processing/OTU.custom.blast.table

# Use Robin's script to calculate the full length pident of each blast result.
cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/taxonomy_processing
Rscript calc_full_length_pident.R OTU.custom.blast.table OTU.custom.blast.table.modified

#############################
# Pull out sequences to run with Silva vs. FreshTrain
#############################
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.above.98 98 TRUE
Rscript filter_seqIDs_by_pident.R OTU.custom.blast.table.modified ids.below.98 98 FALSE
python find_seqIDs_blast_removed.py ../fasta_seqs/16S_seqs_derep.fna OTU.custom.blast.table.modified ids.missing

cat ids.below.98 ids.missing > ids.below.98.all
python create_fastas_given_seqIDs.py ids.above.98 ../fasta_seqs/16S_seqs_derep.fna otus.above.98.fasta
python create_fastas_given_seqIDs.py ids.below.98.all ../fasta_seqs/16S_seqs_derep.fna otus.below.98.fasta

#############################
# Classify sequences
#############################
mothur/mothur "#classify.seqs(fasta=otus.above.98.fasta, template=/home/GLBRCORG/bpeterson26/references/16S/FW.fasta, taxonomy=/home/GLBRCORG/bpeterson26/references/16S/FW.taxonomy, method=wang, probs=T, processors=18)"
mothur/mothur "#classify.seqs(fasta=otus.below.98.fasta, template=/home/GLBRCORG/bpeterson26/references/16S/silva.fasta, taxonomy=/home/GLBRCORG/bpeterson26/references/16S/silva.taxonomy, method=wang, probs=T, processors=18)"

# Clean files
cat otus.above.98.FW.wang.taxonomy otus.below.98.silva.wang.taxonomy > 16S_seqs_derep.taxonomy
sed 's/[[:blank:]]/\;/' <16S_seqs_derep.taxonomy > 16S_seqs_derep.taxonomy.reformatted
mv -f 16S_seqs_derep.taxonomy.reformatted ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/taxonomy/16S_seqs_derep.taxonomy

cd ~/Everglades/dataEdited/16S_from_MG_2018/emirge_processed_data/taxonomy/
sed "s/;/\t/g" 16S_seqs_derep.taxonomy | \
    sed 's/([0-9]*)//g' \
    > 16S_seqs_clean.taxonomy
