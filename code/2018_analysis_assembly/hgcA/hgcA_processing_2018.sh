#!/bin/sh

######################
# code/assembly_analysis/hgcA_processing.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the hgcA sequences
# from our assembly and process them.
######################


############################################
############################################
# hgcA identification
############################################
############################################


######################
# Identify potential hgcA seqs
######################

screen -S EG_hgcA_search

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse/
ORFs=~/Everglades/dataEdited/assemblies/ORFs

mkdir ~/Everglades/dataEdited/2018_analysis_assembly/hgcA
mkdir ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/identification
cd ~/Everglades/dataEdited/2018_analysis_assembly

cat ~/Everglades/metadata/lists/2018_analysis_assembly_list.txt | while read assembly
do
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e hgcA/identification/$assembly\_hgcA_report.txt ]; then
      echo "Search for hgcA in" $assembly
      hmmsearch --tblout hgcA/identification/$assembly\_hgcA.out \
                --cpu 4 \
                --cut_tc \
                ~/references/hgcA/hgcA.hmm \
                $ORFs/$assembly.faa \
                > hgcA/identification/$assembly\_hgcA_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              hgcA/identification/$assembly\_hgcA.out \
              $ORFs/$assembly.faa \
              hgcA/identification/$assembly\_hgcA.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

# exit

# Concatenate all hits and generate list of names
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/identification

cat *_hgcA.faa > hgcA_raw.faa

grep '>' hgcA_raw.faa | \
  sed 's/>//' \
  > hgcA_raw.txt




######################
# Align and quality filter results
######################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/identification

# Align sequences to HMM
hmmalign -o hgcA_raw.sto \
            ~/references/hgcA/hgcA.hmm \
            hgcA_raw.faa

# Convert alignment to fasta format
$scripts/convert_stockhold_to_fasta.py hgcA_raw.sto

# Download hgcA_raw.faa and check it out in Geneious.
# Upload a trimmed hgcA_good.afa fasta file to GLBRC.

grep '>' hgcA_good.afa | \
    sed 's/>//' \
    > hgcA_good.txt

sed 's/-//g' hgcA_good.afa > hgcA_good.faa






############################################
############################################
# Extract gene neighborhood data for hgcA+ scaffolds
############################################
############################################

######################
# First pull out hgcA+ scaffolds
######################

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA
# rm -r scaffolds
mkdir scaffolds


# Pull out scaffold FNA files

cat identification/hgcA_raw.txt \
    | awk -F '_' '{ print $1"_"$2 }' \
    > scaffolds/hgcA_raw_scaffold_list.txt

rm -f scaffolds/hgcA_scaffolds.fna
cat scaffolds/hgcA_raw_scaffold_list.txt | while read scaffold
do

  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "from" $assemblyName

  grep -A 1 $scaffold\$ ~/Everglades/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> scaffolds/hgcA_scaffolds.fna

done


# Pull out GFF entries

rm -f scaffolds/hgcA_scaffolds.gff
cat scaffolds/hgcA_raw_scaffold_list.txt | while read scaffold
do

  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')

  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/Everglades/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> scaffolds/hgcA_scaffolds.gff

done




######################
# Isolate gene neighborhoods
######################

screen -S EG_hgcA_gene_neighborhood

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''

IFS=$'\n'

scripts=~/Everglades/code/generalUse

for hgcA_id in $(cat identification/hgcA_raw.txt)
do

  echo "Working on" $hgcA_id
  scaffold_id=$(echo $hgcA_id | cut -d '_' -f 1-2)
  awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/hgcA_scaffolds.gff > scaffolds/temp_scaffolds.gff

  if [ $(echo $scaffold_id | cut -c5-8 ) == "Mega" ]; then

    echo "MegaHit assembly for" $scaffold_id
    scaffold_num=$(awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/temp_scaffolds.gff | \
                      head -n 1 | \
                      awk -F '\t' '{ print $9 }' | \
                      cut -d"_" -f1 | \
                      sed 's/ID=//' )
    gene_num=$(echo $hgcA_id | cut -d"_" -f3)
    gene_id=$(echo $scaffold_num"_"$gene_num)

  elif [ $(echo $scaffold_id | cut -c5-8 ) == "Meta" ]; then
    echo "metaSPADes assembly for" $scaffold_id
    gene_id=$(echo $hgcA_id | \
                cut -d '_' -f 2-3 | \
                sed 's/^0*//g')

  fi

  echo "Searching for" $gene_id
  python $scripts/gene_neighborhood_extraction.py scaffolds/temp_scaffolds.gff \
                                                  scaffolds/hgcA_scaffolds.fna \
                                                  $gene_id \
                                                  5000 \
                                                  scaffolds/temp_$gene_id

  rm -f scaffolds/temp_scaffolds.gff

done

cd scaffolds
rm -f hgcA_geneNeighborhood.gff hgcA_geneNeighborhood.fna
cat temp_*.gff > hgcA_geneNeighborhood_raw.gff
cat temp_*.fna > hgcA_geneNeighborhood_raw.fna
rm -f *_neighborhood.*



######################
# Keep only samples that passed trimming
######################

IFS=$'\n'

scripts=~/Everglades/code/generalUse
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/

rm -f scaffolds/hgcA_geneNeighborhood.fna
rm -f scaffolds/hgcA_geneNeighborhood.gff

for hgcA_id in $(cat identification/hgcA_good.txt)
do

  echo "Pulling out" $scaffold "fna seqs"
  scaffold=$(echo $hgcA_id | cut -d"_" -f1,2)
  grep -A 1 $scaffold\$ scaffolds/hgcA_geneNeighborhood_raw.fna \
      >> scaffolds/hgcA_geneNeighborhood.fna

  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' scaffolds/hgcA_geneNeighborhood_raw.gff \
      >> scaffolds/hgcA_geneNeighborhood.gff

done



############################################
############################################
# Search for downstream hgcB sequences
############################################
############################################

#########################
# Extract names of downstream seqs
#########################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''

IFS=$'\n'

scripts=~/Everglades/code/generalUse
ORFs=~/Everglades/dataEdited/assemblies/ORFs

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/
mkdir hgcB

rm -f hgcB/downstream_gene_list.txt
python $scripts/retrieve_downstream_gene_name.py \
          identification/hgcA_good.txt \
          scaffolds/hgcA_scaffolds.gff \
          hgcB/downstream_gene_list.txt



#########################
# Extract downstream amino acid sequences
#########################

rm -f hgcB/downstream_genes.faa
cat hgcB/downstream_gene_list.txt | while read gene
do
  assemblyName=$(echo $gene | cut -d "_" -f1)
  echo "Pulling out" $gene "faa entries"
  grep -A 1 $gene$ $ORFs/$assemblyName.faa >> hgcB/downstream_genes.faa
done



#########################
# Search adjacent genes with hgcB HMM
#########################

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/
hmmsearch --tblout hgcB/hgcB.out \
          --cpu 4 \
          -T 30 \
          ~/references/hgcA/hgcB_5M.HMM \
          hgcB/downstream_genes.faa \
          > hgcB/hgcB_report.txt

rm -f hgcB/hgcB.faa
grep -v "#" hgcB/hgcB.out | awk '{ print $1 }' | while read geneID
do
  grep -A 1 $geneID$ hgcB/downstream_genes.faa >> hgcB/hgcB.faa
done

# Align sequences to HMM
hmmalign -o hgcB/hgcB.sto \
            ~/references/hgcA/hgcB_5M.HMM \
            hgcB/hgcB.faa

# Convert alignment to fasta format
$scripts/convert_stockhold_to_fasta.py hgcB/hgcB.sto

# Check sequences in Geneious.
# If they all check out, then just use the hgcB.faa
# file as our final set of hgcB seqs.
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/
grep '>' hgcB/hgcB.faa | \
  sed 's/>//' > hgcB/hgcB.txt


############################################
############################################
# Pull out depth of hgcA+ scaffolds
############################################
############################################

screen -S EG_hgcA_depth

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA
mkdir depth

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  rm -f depth/$metagenome\_hgcA_depth_raw.tsv

  conda activate bioinformatics
  PERL5LIB=""
  cat identification/hgcA_raw.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1 }')
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold ~/Everglades/dataEdited/2018_analysis_assembly/mapping/$metagenome\_to_$assembly.bam \
        >> depth/$metagenome\_hgcA_depth_raw.tsv
  done
  conda deactivate

  echo "Aggregating hgcA depth information for" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python ~/Everglades/code/generalUse/calculate_depth_contigs.py \
            depth/$metagenome\_hgcA_depth_raw.tsv \
            150 \
            depth/$metagenome\_hgcA_depth.tsv
  conda deactivate

  rm -f depth/$metagenome\_hgcA_depth_raw.tsv

done






##################################################
##################################################
# Clustering of sequences
##################################################
##################################################

# Dereplicate hgcA sequences
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA
mkdir dereplication

cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/hgcA_good.faa \
              -o dereplication/hgcA_derep.faa \
              -c 0.97 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/hgcA_derep.faa.clstr > dereplication/hgcA_cluster_faa.tsv























############################################
############################################
# Make phylogenetic tree of sequences
############################################
############################################


#########################
# Find additional reference sequences
#########################

####################
# Retrieve additional reference sequences
####################

screen -S EG_hgcA_references

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''

mkdir ~/Everglades/references/hgcA
mkdir ~/Everglades/references/hgcA/nr_2019_analysis
cd ~/Everglades/references/hgcA/nr_2019_analysis

blastp -query ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/identification/hgcA_good.faa \
        -db /opt/bifxapps/ncbi-blastdb/nr \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 3 \
        -num_threads 3 \
        -out hgcA_ref_blast.txt

# Generate a fasta file of identified sequences
awk -F '\t' '{ print ">"$1"\n"$2 }' hgcA_ref_blast.txt > hgcA_ref_blast.faa
sed 's/gi|[0-9]*|[a-z]*|//' hgcA_ref_blast.faa | sed 's/|//' | sed 's/-//g' > hgcA_ref_blast_edited.faa
mv -f hgcA_ref_blast_edited.faa hgcA_ref_blast.faa

# List of all unique identified sequences
grep '>' hgcA_ref_blast.faa | sed 's/>//' \
                    | sort | uniq \
                    > hgcA_ref_list.txt

# Keep only unique sequence

IFS=$'\n'
for hgcA in $(cat hgcA_ref_list.txt)
do
  grep -A 1 -m 1 $hgcA hgcA_ref_blast.faa >> hgcA_ref_blast_derep.faa
done
mv -f hgcA_ref_blast_derep.faa hgcA_ref_blast.faa

# https://www.ncbi.nlm.nih.gov/books/NBK179288/
# Download the list (hgcA_ref_list.txt) so we can
# get the taxonomic information.
mkdir /Users/benjaminpeterson/Documents/research/Everglades/references/hgcA/nr_2019_analysis
cd /Users/benjaminpeterson/Documents/research/Everglades/references/hgcA/nr_2019_analysis

epost -db protein -input hgcA_ref_list.txt \
    | esummary \
    | xtract -pattern DocumentSummary -element AccessionVersion,Organism > hgcA_ref_taxonomy.tsv



#########################
# Dereplicate reference sequence sets
#########################

cd ~/Everglades/references/hgcA/
cdhit=~/programs/cdhit-master
$cdhit/cd-hit-2d -i hgcA_references_confirmed_methylators.faa \
                  -i2 nr_2019_analysis/hgcA_ref_blast.faa \
                  -o nr_2019_analysis/hgcA_ref_blast_uniq.faa \
                  -c 0.97 \
                  -n 5

cat nr_2019_analysis/hgcA_ref_blast_uniq.faa \
    hgcA_references_confirmed_methylators.faa \
    > nr_2019_analysis/hgcA_ref_for_tree.faa




#########################
# Generate alignment
#########################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA
mkdir phylogeny

# Combine all the sequences we need
cat identification/hgcA_good.faa \
    ~/Everglades/references/hgcA/nr_2019_analysis/hgcA_ref_for_tree.faa \
    ~/5M/dataEdited/binAnalysis/hgcA/identification/hgcA_bins.faa \
    > phylogeny/hgcA_for_phylogeny_raw.faa

# Align sequences using muscle
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/phylogeny
muscle -in hgcA_for_phylogeny_raw.faa \
        -out hgcA_for_phylogeny_raw.afa




#########################
# Generate rough tree
#########################

FastTree hgcA_for_phylogeny_raw.afa \
    > rough_hgcA.tree

# Download this to my computer:
# dataEdited/metagenomes/2018_analysis_assembly/hgcA/phylogeny
# and check it out in Geneious.


#########################
# Generate tree in RAxML
#########################

screen -S EG_hgcA_tree

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/phylogeny

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_masked.afa \
        -n hgcA


























############################################
############################################
# Pull out section of assembly graph for each sequence
############################################
############################################

screen -S assembly_graph
export PATH=/home/GLBRCORG/bpeterson26/miniconda3/bin:$PATH
source activate bandage

mkdir ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/assembly_graphs
cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA

IFS=$'\n'

for geneID in $(cat identification/hgcA_raw.txt)
do

  echo "Working on" $geneID
  assembly=$(echo $geneID | awk -F '_' '{ print $1 }')
  scaffoldName=$(echo $geneID | \
              awk -F '_' '{ print $1"_"$2 }')
  scaffoldID=$(awk -F '\t' -v scaffoldName="$scaffoldName" '$1 == scaffoldName { print $2 }' ~/Everglades/dataEdited/assemblies/scaffolds/renaming_reports/$assembly\_report_file.txt)

  assemblyGraphLocation=~/Everglades/dataEdited/assemblies/assembly_files/$assembly

  grep $scaffoldID $assemblyGraphLocation/assembly_graph_with_scaffolds.gfa | \
    awk '{ print $3 }' | \
    sed 's/-//g' | \
    sed 's/+//g' | \
    sed 's/,/\
/g' > assembly_graphs/nodes_of_interest.txt


  command="Bandage reduce $assemblyGraphLocation/assembly_graph.fastg assembly_graphs/assembly_graph_$scaffoldName.fastg --scope aroundnodes --distance 4 --nodes "

  IFS=$'\n'

  for scaffold_ID in `cat assembly_graphs/nodes_of_interest.txt`
  do
    echo "Adding" $scaffold_ID
    command=$(echo $command$scaffold_ID",")
  done


  command=$(echo $command | sed 's/,$//')
  echo "Running command:"
  echo $command

  bash -c $command

done



######################
# Check out node 1857184499
######################

cd ~/Everglades/dataEdited/2018_analysis_assembly/hgcA/

mkdir reassembly

grep '^P' ~/Everglades/dataEdited/assemblies/assembly_files/KMBP004F/assembly_graph_with_scaffolds.gfa | \
  grep '1857184499' | \
  awk -F '\t' '{ print $2 }'
# scaffold is NODE_26726_length_13659_cov_9.070130
grep 'NODE_26726_length_13659_cov_9.070130' ~/Everglades/dataEdited/assemblies/scaffolds/renaming_reports/KMBP004F_report_file.txt
# Renamed to KMBP004F_000000026726

scaffold=KMBP004F_000000026726
grep -A 1 $scaffold\$ ~/Everglades/dataEdited/assemblies/scaffolds/KMBP004F_assembly.fna \
    > reassembly/node_1857184499.fna
awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/Everglades/dataEdited/assemblies/ORFs/KMBP004F.gff \
    > reassembly/node_1857184499.gff



######################
# Pull out BAM files for each of these
######################

cd ~/Everglades/dataEdited/2018_analysis_assembly/

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

echo -e "KMBP004F_000000026726\nKMBP004F_000000437519\nKMBP004F_000000037185" > hgcA/reassembly/scaffolds_of_interest.txt

IFS=$'\n'

for scaffold in $(cat hgcA/reassembly/scaffolds_of_interest.txt)
do

  echo "Working on" $scaffold

  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')

  samtools view -F 4 \
            mapping/KMBP004F_to_$assemblyName.bam \
            $scaffold \
            > hgcA/reassembly/$scaffold\_reads.sam
  wc -l hgcA/reassembly/$scaffold\_reads.sam

  samtools view -F 8 -f 4 \
            mapping/KMBP004F_to_$assemblyName.bam \
            $scaffold \
            >> hgcA/reassembly/$scaffold\_reads.sam
  wc -l hgcA/reassembly/$scaffold\_reads.sam

done
