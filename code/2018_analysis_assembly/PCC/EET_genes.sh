#!/bin/sh

# code/2018_analysis_assembly/PCC/EET_genes.sh
# Benjamin D. Peterson

##################################################
##################################################
# Identify potential PCCs
##################################################
##################################################

#########################
# Search for MHCs
#########################

screen -S EG_MHCs

cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins
mkdir MHCs
cd MHCs

scripts=~/Everglades/code/generalUse
ORFs=~/Everglades/dataEdited/assemblies/ORFs
assembly_list=~/Everglades/metadata/lists/2018_analysis_assembly_list.txt


# Let's look for proteins with at least 3 haem-binding sites
cat $assembly_list | while read assembly
do
  echo "Searching for MHCs in" $assembly
  $scripts/Find_multiheme_protein.py $ORFs/$assembly.faa 3

  # Move the scripts to the correct directory
  mv $ORFs/$assembly\_3_heme* ./

done

cat *3_heme_list* > MHC_list.txt



#########################
# Extract proteins to either side
#########################

screen -S EG_PCC
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins
mkdir PCC
ORFs=~/Everglades/dataEdited/assemblies/ORFs
IFS=$'\n'
rm -f PCC/adjacent_genes.faa

for gene in $(cat MHCs/MHC_list.txt)
do
  echo "Working on" $gene

  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)

  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)

  preceedingORF=$(echo $scaffold"_"$preceedingORFnumber'$')
  followingORF=$(echo $scaffold"_"$followingORFnumber'$')

  grep -A 1 $preceedingORF $ORFs/$assembly.faa >> PCC/adjacent_genes.faa
  grep -A 1 $followingORF $ORFs/$assembly.faa >> PCC/adjacent_genes.faa

done

# Dereplicate sequence names
cd PCC
grep '>' adjacent_genes.faa | \
  sed 's/>//' | \
  sort | \
  uniq \
  > adjacent_genes_derep_list.txt

# Pull out dereplicated set
rm -f adjacent_genes_derep.faa
cat adjacent_genes_derep_list.txt | while read gene
do
  echo "Working on pulling out" $gene
  grep -A 1 -m 1 $gene adjacent_genes.faa >> adjacent_genes_derep.faa
done










##################################################
##################################################
# Search adjacent proteins for BB-OMPs
##################################################
##################################################

screen -S EG_MHCs

cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/PCC
mkdir BBOMP_HMM

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""


#########################
# Get putative sequences
#########################

# Run the PCC HMM
pcc_omp_HMM=~/Everglades/references/metabolicProteins/BBOMP/pcc_omp.HMM
hmmsearch --tblout BBOMP_HMM/pcc_omp.out \
          -T 50 \
          $pcc_omp_HMM \
          adjacent_genes_derep.faa \
          > BBOMP_HMM/pcc_omp_output.txt
# Pull out the sequences of interest
scripts=~/Everglades/code/generalUse/
python $scripts/extract_protein_hitting_HMM.py \
        BBOMP_HMM/pcc_omp.out \
        adjacent_genes_derep.faa \
        BBOMP_HMM/pcc_omp.faa
# Dereplicate sequences
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/PCC/BBOMP_HMM
cd-hit -i pcc_omp.faa \
        -o pcc_omp_derep.faa \
        -g 1 \
        -c 0.97 \
        -n 5 \
        -d 0

# Align sequences to HMM
hmmalign -o pcc_omp_derep.sto \
            $pcc_omp_HMM \
            pcc_omp_derep.faa
# Convert alignment to fasta format
scripts=~/Everglades/code/generalUse
$scripts/convert_stockhold_to_fasta.py pcc_omp_derep.sto
python $scripts/cleanFASTA.py pcc_omp_derep.afa
mv -f pcc_omp_derep.afa_temp.fasta pcc_omp_derep.afa







##################################################
##################################################
# Generate phylogeny of putative BBOMP proteins
##################################################
##################################################

screen -S EG_BBOMP_phylogeny

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

references=~/Everglades/references/metabolicProteins/BBOMP/pcc_omp.faa

cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/PCC
mkdir BBOMP_phylogeny

#########################
# Generate alignment for phylogeny
#########################
cat BBOMP_HMM/pcc_omp_derep.faa \
    $references \
    > BBOMP_phylogeny/pcc_omp.faa

cd BBOMP_phylogeny
muscle -in pcc_omp.faa \
        -out pcc_omp.afa


#########################
# Generate tree
#########################
FastTree pcc_omp_masked.afa > pcc_omp.tree





##################################################
##################################################
# Extract depths of potential MHCs
##################################################
##################################################

screen -S EG_PCC
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
cd ~/Everglades/dataEdited/2018_analysis_assembly/

# Make list of scaffolds
grep '>' metabolicProteins/PCC/BBOMP_HMM/pcc_omp_derep.afa | \
  sed 's/>//' | \
  cut -d"_" -f1,2 \
  > metabolicProteins/PCC/BBOMP_HMM/pcc_omp_derep_scaffolds.txt

# Pull out all depths
cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  rm -f metabolicProteins/PCC/BBOMP_HMM/$metagenome\_depth_raw.tsv
  conda activate bioinformatics
  PERL5LIB=""

  # Calculate depth over each residue in each scaffold
  cat metabolicProteins/PCC/BBOMP_HMM/pcc_omp_derep_scaffolds.txt | while read scaffold
  do
    assembly=$(echo $scaffold | awk -F '_' '{ print $1 }')
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold mapping/$metagenome\_to_$assembly.bam \
        >> metabolicProteins/PCC/BBOMP_HMM/$metagenome\_depth_raw.tsv
  done
  conda deactivate

  # Average the depth of each residue over the entire scaffold
  echo "Aggregating" $scaffold "depth information for" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python ~/Everglades/code/generalUse/calculate_depth_contigs.py \
            metabolicProteins/PCC/BBOMP_HMM/$metagenome\_depth_raw.tsv \
            150 \
            metabolicProteins/PCC/BBOMP_HMM/$metagenome\_depth.tsv
  conda deactivate

done

rm -f metabolicProteins/PCC/BBOMP_HMM/*_depth_raw.tsv
