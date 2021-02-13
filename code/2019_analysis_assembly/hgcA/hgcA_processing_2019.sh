#!/bin/sh

######################
# code/assembly_analysis/hgcA_processing_2019.sh
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
cd ~/Everglades/metadata/lists/


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

mkdir ~/Everglades/dataEdited/2019_analysis_assembly/hgcA
mkdir ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/identification
cd ~/Everglades/dataEdited/2019_analysis_assembly


cat ~/Everglades/metadata/lists/2019_analysis_assembly_list.txt | while read assembly
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
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/identification

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

cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/identification

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

cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA
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
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA
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

  if [ $(echo $scaffold_id | cut -c7-10 ) == "Mega" ]; then

    echo "MegaHit assembly for" $scaffold_id
    scaffold_num=$(awk -F '\t' -v scaffold_id="$scaffold_id" '$1 == scaffold_id { print $0 }' scaffolds/temp_scaffolds.gff | \
                      head -n 1 | \
                      awk -F '\t' '{ print $9 }' | \
                      cut -d"_" -f1 | \
                      sed 's/ID=//' )
    gene_num=$(echo $hgcA_id | cut -d"_" -f3)
    gene_id=$(echo $scaffold_num"_"$gene_num)

  elif [ $(echo $scaffold_id | cut -c7-10 ) == "Meta" ]; then
    echo "metaSPADes assembly for" $scaffold_id
    gene_id=$(echo $hgcA_id | \
                cut -d '_' -f 2-3 | \
                sed 's/^0*//g')

  else
    echo "You fucked up" $scaffold_id "isn't ID'd as MegaHit or metaSPADEs."
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
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/
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
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/
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

cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

scripts=~/Everglades/code/generalUse
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/
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
# I ended up removing a bunch of them.
# Cleaned version is hgcB_good.afa.
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/
sed 's/-//g' hgcB/hgcB_good.afa > hgcB/hgcB_good.faa

grep '>' hgcB/hgcB_good.faa | \
  sed 's/>//' > hgcB/hgcB_good.txt


############################################
############################################
# Pull out depth of hgcA+ scaffolds
############################################
############################################

screen -S EG_hgcA_depth

cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA
mkdir depth

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cat ~/Everglades/metadata/lists/2019_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  rm -f depth/$metagenome\_hgcA_depth_raw.tsv

  conda activate bioinformatics
  PERL5LIB=""
  cat identification/hgcA_raw.txt | while read gene
  do
    scaffold=$(echo $gene | awk -F '_' '{ print $1"_"$2 }')
    assembly=$(echo $gene | awk -F '_' '{ print $1 }')
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold ~/Everglades/dataEdited/2019_analysis_assembly/mapping/$metagenome\_to_$assembly.bam \
        >> depth/$metagenome\_hgcA_depth_raw.tsv
  done
  conda deactivate

  echo "Aggregating hgcA depth information for" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python ~/Everglades/code/generalUse/calculate_depth_length_contigs.py \
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
cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA
mkdir dereplication


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/hgcA_good.faa \
              -o dereplication/hgcA_derep_divergent.faa \
              -c 0.80 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/hgcA_derep_divergent.faa.clstr > dereplication/hgcA_divergent_cluster_faa.tsv
# Download this.


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/hgcA_good.faa \
              -o dereplication/hgcA_derep_96.faa \
              -c 0.96 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/hgcA_derep_96.faa.clstr > dereplication/hgcA_cluster_faa_96.tsv
# Download this.


cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i identification/hgcA_good.faa \
              -o dereplication/hgcA_derep_94.faa \
              -c 0.94 \
              -n 5 \
              -d 0
$cdhit/clstr2txt.pl dereplication/hgcA_derep_94.faa.clstr > dereplication/hgcA_cluster_faa_94.tsv
# Download this.



# Final set of hgcA sequences:
grep -A 1 -f dereplication/hgcA_derep_list.txt identification/hgcA_good.faa \
  > hgcA.faa
grep -v "\-\-" hgcA.faa > hgcA_clean.faa
mv -f hgcA_clean.faa hgcA.faa
cp -f dereplication/hgcA_derep_list.txt ./hgcA.txt


cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/
echo -e "hgcA\tclusterRep" > derep_key.tsv
cat identification/hgcA_good.txt | while read hgcA
do
  clusterID=`awk -F '\t' -v hgcA="$hgcA" ' $1 == hgcA { print $2 }' dereplication/hgcA_cluster_faa_94.tsv`
  clusterRep=`awk -F '\t' -v clusterID="$clusterID" ' ($2 == clusterID) && ($5 == 1) { print $1 }' dereplication/hgcA_cluster_faa_94.tsv`
  echo -e $hgcA'\t'$clusterRep >> derep_key.tsv
done





############################################
############################################
# Classify hgcA seqs with pplacer workflow
############################################
############################################

screen -S EG_hgcA_pplacer
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
references=~/Everglades/references/hgcA
workingDirectory=~/Everglades/dataEdited/2019_analysis_assembly/hgcA
scripts=~/Everglades/code/generalUse
mkdir $workingDirectory/classification

# Generate fasta file of reference alignment
cd $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg
$scripts/convert_stockhold_to_fasta.py Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.stockholm
mv Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afackholm $workingDirectory/classification/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa

# Generate alignment of sequences of interest
cd $workingDirectory
muscle -in hgcA.faa \
        -out classification/hgcA_muscle.afa
cd classification/

# Combine the alignments of seqs from this study and references
muscle -profile -in1 hgcA_muscle.afa \
        -in2 Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.afa \
        -out hgcA_for_classification.afa
# Convert to stockholm format
python $scripts/convert_fasta_to_stockholm.py hgcA_for_classification.afa
conda deactivate

# Run pplacer
conda activate hgcA_classifier
pplacer -p \
        --keep-at-most 1 \
        --max-pend 1 \
        -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/ \
        hgcA_for_classification.sto

# Make sqlite database of reference
rppr prep_db \
      --sqlite Hg_MATE_classify \
      -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg

# Generate taxonomic assignments using guppy
guppy classify -c $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg \
               --pp \
               --sqlite Hg_MATE_classify \
               hgcA_for_classification.jplace

# Save out this data to csv
guppy to_csv --point-mass \
              --pp \
              -o hgcA_classification.csv \
              hgcA_for_classification.jplace

# Visualize placements on FastTree
guppy tog --pp \
          -o hgcA_classification.nwk \
          hgcA_for_classification.jplace




############################################
############################################
# Make phylogenetic tree
############################################
############################################

####################
# Make phylogenetic tree with Hg-MATE database for references
####################

screen -S EG_hgcA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''

# Generate alignment
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
workingDirectory=~/Everglades/dataEdited/2019_analysis_assembly/hgcA/phylogeny
scripts=~/Everglades/code/generalUse/
mkdir $workingDirectory
cd $workingDirectory

muscle -profile \
        -in1 ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/classification/hgcA_muscle.afa \
        -in2 $references/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full.refpkg/Hg-MATE-Db.v1.ISOCELMAG_HgcA_full_align-bmge.fasta \
        -out hgcA_for_phylogeny_raw.afa
# Generate rough tree
FastTree hgcA_for_phylogeny_raw.afa \
    > rough_hgcA.tree
# Download this to my local computer.

# Upload the list of sequences to remove
python $scripts/remove_fasta_seqs_using_list_of_headers.py hgcA_for_phylogeny_raw.afa \
                                                            seqs_to_remove.txt \
                                                            hgcA_for_phylogeny.afa


#########################
# Generate tree in RAxML
#########################

screen -S EG_hgcA_tree

cd ~/Everglades/dataEdited/2019_analysis_assembly/hgcA/phylogeny

raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
#$raxml -f a \
#        -p 283976 \
#        -m PROTGAMMAAUTO \
#        -N autoMRE \
#        -x 2381 \
#        -T 20 \
#        -s hgcA_for_phylogeny_masked.afa \
#        -n hgcA
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s hgcA_for_phylogeny_masked.afa.reduced \
        -n hgcA
