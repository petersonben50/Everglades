#!/bin/sh

######################
# code/assembly_analysis/merB/merB_processing_2019.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to pull out the merB sequences
# from our assembly and process them.
######################


############################################
############################################
# merB identification
############################################
############################################
cd ~/Everglades/metadata/lists/


######################
# Identify potential merB seqs
######################

screen -S EG_merB_search

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse/
ORFs=~/Everglades/dataEdited/assemblies/ORFs

mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB
mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB/identification
cd ~/Everglades/dataEdited/2019_analysis_assembly


cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e merB/identification/$assembly\_merB_report.txt ]; then
      echo "Search for merB in" $assembly
      hmmsearch --tblout merB/identification/$assembly\_merB.out \
                --cpu 4 \
                -E 0.00001 \
                ~/references/merB/merB_christakis.hmm \
                $ORFs/$assembly.faa \
                > merB/identification/$assembly\_merB_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              merB/identification/$assembly\_merB.out \
              $ORFs/$assembly.faa \
              merB/identification/$assembly\_merB.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

# Concatenate all results
cd merB/identification
head -n 2 Sed992Mega19_merB.out | \
  tail -n 1 | \
  > all_merB_outputs.txt
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  grep -v '^#' $assembly\_merB.out >> all_merB_outputs.txt
done
# exit





######################
# Identify potential merB seqs
######################

screen -S EG_merB_pfam_search

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse/
ORFs=~/Everglades/dataEdited/assemblies/ORFs

mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB_pfam
mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB_pfam/identification
cd ~/Everglades/dataEdited/2019_analysis_assembly


cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e merB_pfam/identification/$assembly\_merB_pfam_report.txt ]; then
      echo "Search for merB_pfam in" $assembly
      hmmsearch --tblout merB_pfam/identification/$assembly\_merB_pfam.out \
                --cpu 4 \
                -E 0.00001 \
                ~/references/merB/MerB.hmm \
                $ORFs/$assembly.faa \
                > merB_pfam/identification/$assembly\_merB_pfam_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              merB_pfam/identification/$assembly\_merB_pfam.out \
              $ORFs/$assembly.faa \
              merB_pfam/identification/$assembly\_merB_pfam.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

# Concatenate all results
cd merB_pfam/identification
head -n 2 Sed992Mega19_merB_pfam.out | \
  tail -n 1 | \
  > all_merB_pfam_outputs.txt
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  grep -v '^#' $assembly\_merB_pfam.out >> all_merB_pfam_outputs.txt
done
# exit


# Concatenate all hits and generate list of names
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/identification
cat *_merB.faa > merB_raw.faa
grep '>' merB_raw.faa | \
  sed 's/>//' \
  > merB_raw.txt



  ######################
  # Align and quality filter results from Christakis paper
  ######################

  source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
  conda activate bioinformatics
  PYTHONPATH=''
  PERL5LIB=''

  scripts=~/Everglades/code/generalUse

  cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/identification
  cat ~/references/merB/merB_refs.faa merB_raw.faa > merB_raw_with_refs.faa

  # Align sequences to HMM
  hmmalign -o merB_raw.sto \
              ~/references/merB/merB_christakis.hmm \
              merB_raw_with_refs.faa
  muscle -in merB_raw_with_refs.faa \
        -out merB_raw_with_refs_muscle.afa
  # Convert alignment to fasta format
  $scripts/convert_stockhold_to_fasta.py merB_raw.sto

  # Download merB_raw.faa and check it out in Geneious.
  # Upload a trimmed merB_good.afa fasta file to GLBRC.

  grep '>' merB_good.afa | \
      sed 's/>//' \
      > merB_good.txt

  sed 's/-//g' merB_good.afa > merB_good.faa



  cd ~/Documents/research/Everglades/dataEdited/assembly_analysis/merB/identification/outputs
  ls merB*afa | sed 's/.afa//' | while read cluster
  do
    grep ">" $cluster.afa | sed 's/>//' > $cluster.txt
  done






######################
# Identify potential merB seqs
######################

screen -S EG_merB_final_search

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse/
ORFs=~/Everglades/dataEdited/assemblies/ORFs

mv ~/Everglades/dataEdited/2019_analysis_assembly/merB \
  ~/Everglades/dataEdited/2019_analysis_assembly/merB_testing


mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB
mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB/identification
cd ~/Everglades/dataEdited/2019_analysis_assembly


cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  if [ -e $ORFs/$assembly.faa ]; then
    if [ ! -e merB/identification/$assembly\_merB_report.txt ]; then
      echo "Search for merB in" $assembly
      hmmsearch --tblout merB/identification/$assembly\_merB.out \
                --cpu 4 \
                -T 75 \
                ~/references/merB/merB_christakis.hmm \
                $ORFs/$assembly.faa \
                > merB/identification/$assembly\_merB_report.txt
      python $scripts/extract_protein_hitting_HMM.py \
              merB/identification/$assembly\_merB.out \
              $ORFs/$assembly.faa \
              merB/identification/$assembly\_merB.faa
    else
      echo "Search is already done in" $assembly
    fi
  else
    echo "Genes aren't predicted for" $assembly
  fi
done

# Concatenate all results
cd merB/identification
head -n 2 Sed992Mega19_merB.out | \
  tail -n 1 | \
  > all_merB_outputs.txt
cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
do
  grep -v '^#' $assembly\_merB.out >> all_merB_outputs.txt
done
# exit


# Concatenate all hits and generate list of names
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/identification
cat *_merB.faa > merB_raw.faa
grep '>' merB_raw.faa | \
  sed 's/>//' \
  > merB_raw.txt







############################################
############################################
# Extract gene neighborhood data for merB+ scaffolds
############################################
############################################

######################
# First pull out merB+ scaffolds
######################

cd ~/Everglades/dataEdited/2019_analysis_assembly/merB
# rm -r scaffolds
mkdir scaffolds


# Pull out scaffold FNA files
cat identification/merB_raw.txt \
    | awk -F '_' '{ print $1"_"$2 }' \
    > scaffolds/merB_raw_scaffold_list.txt

rm -f scaffolds/merB_scaffolds.fna
cat scaffolds/merB_raw_scaffold_list.txt | while read scaffold
do

  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')
  echo "Pulling out" $scaffold "from" $assemblyName

  grep -A 1 $scaffold\$ ~/Everglades/dataEdited/assemblies/scaffolds/$assemblyName\_assembly.fna \
      >> scaffolds/merB_scaffolds.fna

done


# Pull out GFF entries

rm -f scaffolds/merB_scaffolds.gff
cat scaffolds/merB_raw_scaffold_list.txt | while read scaffold
do

  assemblyName=$(echo $scaffold | awk -F '_' '{ print $1 }')

  echo "Pulling out" $scaffold "GFF entries"
  awk -v scaffold="$scaffold" '{ if ($1 == scaffold) print }' ~/Everglades/dataEdited/assemblies/ORFs/$assemblyName.gff \
      >> scaffolds/merB_scaffolds.gff

done




######################
# Isolate gene neighborhoods
######################

screen -S EG_merB_gene_neighborhood
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate py_viz
PYTHONPATH=''
IFS=$'\n'
scripts=~/Everglades/code/generalUse

cat identification/merB_raw.txt | while read merB_id
do

  echo "Working on" $merB_id
  scaffold_id=$(echo $merB_id | cut -d '_' -f 1-2)
  awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/merB_scaffolds.gff > scaffolds/temp_scaffolds.gff

  if [ $(echo $scaffold_id | cut -c7-10 ) == "Mega" ]; then

    echo "MegaHit assembly for" $scaffold_id
    scaffold_num=$(awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/temp_scaffolds.gff | \
                      head -n 1 | \
                      awk -F '\t' '{ print $9 }' | \
                      cut -d"_" -f1 | \
                      sed 's/ID=//' )
    gene_num=$(echo $merB_id | cut -d"_" -f3)
    gene_id=$(echo $scaffold_num"_"$gene_num)

  elif [ $(echo $scaffold_id | cut -c7-10 ) == "Meta" ]; then
    echo "metaSPADes assembly for" $scaffold_id
    gene_id=$(echo $merB_id | \
                cut -d '_' -f 2-3 | \
                sed 's/^0*//g')

  else
    echo "You fucked up" $scaffold_id "isn't ID'd as MegaHit or metaSPADEs."
  fi

  echo "Searching for" $gene_id
  python $scripts/gene_neighborhood_extraction.py scaffolds/temp_scaffolds.gff \
                                                  scaffolds/merB_scaffolds.fna \
                                                  $gene_id \
                                                  5000 \
                                                  scaffolds/temp_$gene_id

  rm -f scaffolds/temp_scaffolds.gff
done

cd scaffolds
rm -f merB_geneNeighborhood.gff merB_geneNeighborhood.fna
cat temp_*.gff > merB_geneNeighborhood_raw.gff
cat temp_*.fna > merB_geneNeighborhood_raw.fna
rm -f *_neighborhood.*




############################################
############################################
# Pull out depth of merB+ scaffolds
############################################
############################################

screen -S EG_merB_depth

cd ~/Everglades/dataEdited/2019_analysis_assembly/merB
mkdir depth

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  rm -f depth/$metagenome\_merB_depth_raw.tsv

  conda activate bioinformatics
  PERL5LIB=""
  cat identification/merB_raw.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1 }')
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold ~/Everglades/dataEdited/2019_analysis_assembly/mapping/$metagenome\_to_$assembly.bam \
        >> depth/$metagenome\_merB_depth_raw.tsv
  done
  conda deactivate

  echo "Aggregating merB depth information for" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python ~/Everglades/code/generalUse/calculate_depth_length_contigs.py \
            depth/$metagenome\_merB_depth_raw.tsv \
            150 \
            depth/$metagenome\_merB_depth.tsv
  conda deactivate

  rm -f depth/$metagenome\_merB_depth_raw.tsv

done






##################################################
##################################################
# Clustering of sequences
##################################################
##################################################

# Dereplicate merB sequences
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB
mkdir dereplication


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/merB_good.faa \
              -o dereplication/merB_derep_divergent.faa \
              -c 0.80 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/merB_derep_divergent.faa.clstr > dereplication/merB_divergent_cluster_faa.tsv
# Download this.


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/merB_good.faa \
              -o dereplication/merB_derep_96.faa \
              -c 0.96 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/merB_derep_96.faa.clstr > dereplication/merB_cluster_faa_96.tsv
# Download this.


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/merB_good.faa \
              -o dereplication/merB_derep_94.faa \
              -c 0.94 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/merB_derep_94.faa.clstr > dereplication/merB_cluster_faa_94.tsv
# Download this.



# Final set of merB sequences:
grep -A 1 -f dereplication/merB_derep_list.txt identification/merB_good.faa \
  > merB.faa
grep -v "\-\-" merB.faa > merB_clean.faa
mv -f merB_clean.faa merB.faa
cp -f dereplication/merB_derep_list.txt ./merB.txt


cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/
echo -e "merB\tclusterRep" > derep_key.tsv
cat identification/merB_good.txt | while read merB
do
  clusterID=`awk -F '\t' -v merB="$merB" ' $1 == merB { print $2 }' dereplication/merB_cluster_faa_94.tsv`
  clusterRep=`awk -F '\t' -v clusterID="$clusterID" ' ($2 == clusterID) && ($5 == 1) { print $1 }' dereplication/merB_cluster_faa_94.tsv`
  echo -e $merB'\t'$clusterRep >> derep_key.tsv
done





############################################
############################################
# Classify merB seqs with pplacer workflow
############################################
############################################

screen -S EG_merB_pplacer
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
references=~/Everglades/references/merB
workingDirectory=~/Everglades/dataEdited/2019_analysis_assembly/merB
scripts=~/Everglades/code/generalUse
mkdir $workingDirectory/classification

# Generate fasta file of reference alignment
cd $references/Hg-MATE-Db.v1.ISOCELMAG_merB_full.refpkg
$scripts/convert_stockhold_to_fasta.py Hg-MATE-Db.v1.ISOCELMAG_merB_full.stockholm
mv Hg-MATE-Db.v1.ISOCELMAG_merB_full.afackholm $workingDirectory/classification/Hg-MATE-Db.v1.ISOCELMAG_merB_full.afa

# Generate alignment of sequences of interest
cd $workingDirectory
muscle -in merB.faa \
        -out classification/merB_muscle.afa
cd classification/

# Combine the alignments of seqs from this study and references
muscle -profile -in1 merB_muscle.afa \
        -in2 Hg-MATE-Db.v1.ISOCELMAG_merB_full.afa \
        -out merB_for_classification.afa
# Convert to stockholm format
python $scripts/convert_fasta_to_stockholm.py merB_for_classification.afa
conda deactivate

# Run pplacer
conda activate merB_classifier
pplacer -p \
        --keep-at-most 1 \
        --max-pend 1 \
        -c $references/Hg-MATE-Db.v1.ISOCELMAG_merB_full.refpkg/ \
        merB_for_classification.sto

# Make sqlite database of reference
rppr prep_db \
      --sqlite Hg_MATE_classify \
      -c $references/Hg-MATE-Db.v1.ISOCELMAG_merB_full.refpkg

# Generate taxonomic assignments using guppy
guppy classify -c $references/Hg-MATE-Db.v1.ISOCELMAG_merB_full.refpkg \
               --pp \
               --sqlite Hg_MATE_classify \
               merB_for_classification.jplace

# Save out this data to csv
guppy to_csv --point-mass \
              --pp \
              -o merB_classification.csv \
              merB_for_classification.jplace

# Visualize placements on FastTree
guppy tog --pp \
          -o merB_classification.nwk \
          merB_for_classification.jplace




############################################
############################################
# Make phylogenetic tree
############################################
############################################

####################
# Prepped alignment, preliminary tree
####################

screen -S EG_merB_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''

# Generate alignment
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
workingDirectory=~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny
scripts=~/Everglades/code/generalUse/
mkdir $workingDirectory
cd $workingDirectory

muscle -profile \
        -in1 ~/Everglades/dataEdited/2019_analysis_assembly/merB/classification/merB_muscle.afa \
        -in2 $references/Hg-MATE-Db.v1.ISOCELMAG_merB_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_merB_full_align-bmge.fasta \
        -out merB_for_phylogeny_raw.afa
# Generate rough tree
FastTree merB_for_phylogeny_raw.afa \
    > rough_merB.tree
# Download this to my local computer.

# Upload the list of sequences to remove
python $scripts/remove_fasta_seqs_using_list_of_headers.py merB_for_phylogeny_raw.afa \
                                                            seqs_to_remove.txt \
                                                            merB_for_phylogeny.afa


#########################
# Generate tree in RAxML
#########################

screen -S EG_merB_tree

cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
#$raxml -f a \
#        -p 283976 \
#        -m PROTGAMMAAUTO \
#        -N autoMRE \
#        -x 2381 \
#        -T 20 \
#        -s merB_for_phylogeny_masked.afa \
#        -n merB
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s merB_for_phylogeny_masked.afa.reduced \
        -n merB







#########################
# Generate tree for publications with subset of references using RAxML
#########################

# Done locally
Everglades
cp /Users/benjaminpeterson/Documents/research/Hg_MATE/versions/v1.01142021/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    references
mkdir dataEdited/assembly_analysis/merB/phylogeny/final
mkdir dataEdited/assembly_analysis/merB/phylogeny/final/refs
finalFolderRefs=dataEdited/assembly_analysis/merB/phylogeny/final/refs
rm -f $finalFolderRefs/HgMate_reference_seqs_to_use.faa
cut -d"_" -f1-3 $finalFolderRefs/reference_names_to_use.txt | while read reference_name
do
  grep -A 1 \
    $reference_name \
    references/Hg-MATE-Db.v1.01142021_ISOCELMAG_Hgc.fas \
    >> $finalFolderRefs/HgMate_reference_seqs_to_use.faa
done
grep ">" $finalFolderRefs/HgMate_reference_seqs_to_use.faa | \
  sed 's/>//' | tr -d '[:blank:]' \
  > $finalFolderRefs/IDed_seqs.txt
wc -l $finalFolderRefs/*
cat $finalFolderRefs/IDed_seqs.txt $finalFolderRefs/reference_names_to_use.txt | \
  sort

# Done on GLBRC
screen -S EG_merB_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
mkdir ~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny/final
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny/final

# Get seqs from study
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/classification
python ~/Everglades/code/generalUse/cleanFASTA.py merB_muscle.afa
sed 's/-//g' merB_muscle.afa_temp.fasta > ~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny/final/merB.faa


# Upload needed references
screen -S RAxML_merB
cd ~/Everglades/dataEdited/2019_analysis_assembly/merB/phylogeny/final
cat HgMate_reference_seqs_to_use.faa \
    merB.faa \
    > merB_for_tree_final.faa

# Generate alignment
muscle -in merB_for_tree_final.faa \
        -out merB_for_tree_final.afa

# Use trimal to trim alignment
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
trimal -in merB_for_tree_final.afa \
        -out merB_for_tree_final_masked.afa \
        -gt 0.5
conda deactivate

# Upload masked alignment (50% gaps)
# Then run RAxML to generate tree
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s merB_for_tree_final_masked.afa \
        -n merB
