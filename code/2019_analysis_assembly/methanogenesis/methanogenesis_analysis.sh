#!/bin/sh

###############################
# code/2019_analysis_assembly/methanogenesis/methanogenesis_analysis.sh
# Benjamin D. Peterson
###############################


###############################
###############################
# Initial mcrA tree
###############################
###############################

screen -S mcrA_refs
# Add Daan's references
cd ~/Everglades/references/
mkdir metabolicProteins/methanogenesis
# Add his references here

# Get set up for generating tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins
mkdir methanogenesis
mkdir methanogenesis/mcrA
mkdir methanogenesis/mcrA/phylogeny

###############################
# Pull out references from refseq
###############################
blastp -query derep/mcrA_clean.faa \
        -db ~/references/ncbi_db/refseq/refseq_protein \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out methanogenesis/mcrA/phylogeny/refseq_protein_mcrA.txt

# Generate a fasta file of identified sequences
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/
awk -F '\t' '{ print ">"$1"\n"$2 }' refseq_protein_mcrA.txt | \
  sed 's/ref|//' | \
  sed 's/|//' | \
  sed 's/-//g' > refseq_protein_mcrA.faa

# List of all unique identified sequences
grep '>' refseq_protein_mcrA.faa | sed 's/>//' \
                    | sort | uniq \
                    > refseq_protein_mcrA_list.txt

# Keep only unique sequences
IFS=$'\n'
rm -f refseq_protein_mcrA_derep.faa
for mcrA in $(cat refseq_protein_mcrA_list.txt)
do
  grep -A 1 -m 1 $mcrA refseq_protein_mcrA.faa >> refseq_protein_mcrA_derep.faa
done
mv -f refseq_protein_mcrA_derep.faa refseq_protein_mcrA.faa


# Dereplicate refseq sequences
cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i refseq_protein_mcrA.faa \
              -o refseq_protein_mcrA_derep.faa \
              -c 0.97 \
              -n 5 \
              -d 0
$cdhit/cd-hit-2d -i ~/Everglades/references/metabolicProteins/methanogenesis/mcrA_Fig2_clust90.faa \
                  -i2 refseq_protein_mcrA_derep.faa \
                  -o refseq_protein_mcrA_derep_unique.faa \
                  -c 0.97 \
                  -n 5
grep '>' refseq_protein_mcrA_derep_unique.faa | sed 's/>//' \
                    | sort | uniq \
                    > refseq_protein_mcrA_list_derep_unique.txt


###############################
# Generate mcrA tree
###############################
# Concatenate all needed files.
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/
cat ~/Everglades/references/metabolicProteins/methanogenesis/mcrA_Fig2_clust90.faa \
    derep/mcrA_clean.faa \
    methanogenesis/mcrA/phylogeny/refseq_protein_mcrA_derep_unique.faa \
    > methanogenesis/mcrA/phylogeny/mcrA_for_phylogeny.faa

# Generate alignment for tree
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/
sed -i 's/*//' mcrA_for_phylogeny.faa
muscle -in mcrA_for_phylogeny.faa \
        -out mcrA_for_phylogeny.afa

# Generate tree
FastTree mcrA_for_phylogeny_masked.afa > mcrA_for_phylogeny.tree

# Get reference metadata
Everglades
cd dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny
epost -db protein -input refseq_protein_mcrA_list_derep_unique.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism > ref_mcrA_metadata.tsv
