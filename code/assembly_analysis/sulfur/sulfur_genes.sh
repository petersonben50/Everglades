#!/bin/sh

######################
# code/2019_analysis_assembly/sulfur/sulfur_genes.sh
# Benjamin D. Peterson
######################

############################################
############################################
# Sulfate reducers
############################################
############################################

######################
# dsrA
######################

screen -S EG_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/Everglades/code/generalUse
# Set up directory
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins
mkdir sulfur
mkdir sulfur/dsrA

# Copy over my sequences
cp derep/dsrA.afa sulfur/dsrA/dsrA.afa
# Copy in reference sequences
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa


# Generate ML tree
#screen -S EG_dsrA_tree
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s dsrA_phylogeny_masked.afa \
        -n dsrA
FastTree dsrA_phylogeny_masked.afa \
    > dsrA_phylogeny_masked.tree


############################################
############################################
# Reductive dsrA phylogeny
############################################
############################################

screen -S EG_dsrA_references
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA
mkdir dsrA_reductive_phylogeny
# Add in dsrA_red_list.txt file

rm -f dsrA_reductive_phylogeny/dsrA.afa
cat dsrA_reductive_phylogeny/dsrA_red_list.txt | while read dsrA
do
  grep -A 1 -m 1 $dsrA dsrA.afa >> dsrA_reductive_phylogeny/dsrA.afa
done
cd dsrA_reductive_phylogeny
sed 's/-//g' dsrA.afa > dsrA.faa


######################
# Retrieve dsrA references from refseq_protein
######################

# Search for references in the refseq_protein database
blastp -query dsrA.faa \
        -db ~/references/ncbi_db/refseq/refseq_protein \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out refseq_protein_dsrA.txt


# Generate a fasta file of identified sequences
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/dsrA_reductive_phylogeny
awk -F '\t' '{ print ">"$1"\n"$2 }' refseq_protein_dsrA.txt > refseq_protein_dsrA.faa
# List of all unique identified sequences
grep '>' refseq_protein_dsrA.faa | sed 's/>//' \
                    | sort | uniq \
                    > refseq_protein_dsrA_list.txt

# Keep only unique sequences
IFS=$'\n'
for dsrA in $(cat refseq_protein_dsrA_list.txt)
do
  grep -A 1 -m 1 $dsrA refseq_protein_dsrA.faa >> refseq_protein_dsrA_derep.faa
done
mv -f refseq_protein_dsrA_derep.faa refseq_protein_dsrA.faa

# Clean up
sed -i 's/-//g' refseq_protein_dsrA.faa
sed -i 's/ref|//' refseq_protein_dsrA.faa
sed -i 's/|//' refseq_protein_dsrA.faa

# Dereplicate sequences
cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i refseq_protein_dsrA.faa \
              -o refseq_protein_dsrA_derep.faa \
              -c 0.97 \
              -n 5 \
              -d 0


######################
# Retrieve dsrA references from nr
######################

cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/dsrA_reductive_phylogeny
# Search for references in the nr database
blastp -query dsrA.faa \
        -db ~/references/ncbi_db/nr/nr \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out nr_dsrA.txt
# Generate a fasta file of identified sequences
awk -F '\t' '{ print ">"$1"\n"$2 }' nr_dsrA.txt > nr_dsrA.faa

# List of all unique identified sequences
grep '>' nr_dsrA.faa | sed 's/>//' \
                    | sort | uniq \
                    > nr_dsrA_list.txt

# Keep only unique sequences
IFS=$'\n'
for dsrA in $(cat nr_dsrA_list.txt)
do
  grep -A 1 -m 1 $dsrA nr_dsrA.faa >> nr_dsrA_derep.faa
done
mv -f nr_dsrA_derep.faa nr_dsrA.faa


# Clean up fasta file
sed -i 's/-//g' nr_dsrA.faa
sed -i 's/[a-z]*|//' nr_dsrA.faa
sed -i 's/|//' nr_dsrA.faa

# Dereplicate nr sequences
cdhit=~/programs/cdhit-master
$cdhit/cd-hit -g 1 \
              -i nr_dsrA.faa \
              -o nr_dsrA_derep.faa \
              -c 0.97 \
              -n 5 \
              -d 0

# Dereplicate nr seqs against the refseq set
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/dsrA_reductive_phylogeny
cdhit=~/programs/cdhit-master
$cdhit/cd-hit-2d -i refseq_protein_dsrA_derep.faa \
                  -i2 nr_dsrA_derep.faa \
                  -o nr_dsrA_uniq.faa \
                  -c 0.97 \
                  -n 5

# Concatenate all sequences
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsrA/dsrA_reductive_phylogeny
cat dsrA.faa \
    refseq_protein_dsrA_derep.faa \
    nr_dsrA_uniq.faa \
    > dsrA_for_taxonomy_raw.faa


# Retreive metadata for reference sequences
grep --no-filename '>' nr_dsrA_uniq.faa refseq_protein_dsrA_derep.faa | \
  sed 's/>//' \
  > ref_dsrA_list.txt


# On local computer
Everglades
cd dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/reductive_dsrA_phylogeny/
epost -db protein -input ref_dsrA_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism,TaxId > ref_dsrA_metadata.tsv

awk -F '\t' '{ print $3 }' ref_dsrA_metadata.tsv | sort | uniq > ref_taxID.txt

epost -db taxonomy -input ref_taxID.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element TaxId,Division,Genus,Species |
    awk -F '\t' -v OFS=',' '{ print $1,$2,$3,$4 }' > ref_dsrA_taxInfo.tsv



######################
# Generate tree
######################

# Generate alignment
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
muscle -in dsrA_for_taxonomy_raw.faa \
        -out dsrA_for_taxonomy_raw.afa

# Generate rough tree
FastTree dsrA_for_taxonomy_raw.afa \
    > rough_dsrA.tree

# Mask at 50% gaps

# Generate good tree
screen -S EG_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
cd ~/Everglades/dataEdited/2019_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_reductive_phylogeny
raxmlHPC-PTHREADS -f a \
                  -p 283976 \
                  -m PROTGAMMAAUTO \
                  -N autoMRE \
                  -x 2381 \
                  -T 20 \
                  -s dsrA_for_taxonomy_masked.afa \
                  -n dsrA
