#!/bin/sh

#########################
# code/2018_analysis_assembly/metabolic_proteins.sh
# Benjamin D. Peterson
#########################


#########################
# Identification of metabolic proteins
#########################

screen -S EG_metabolic_HMMs

cd ~/Everglades/dataEdited/2018_analysis_assembly/
mkdir metabolicProteins
mkdir metabolicProteins/identification

scripts=~/Everglades/code/generalUse
metabolic_HMMs=~/Everglades/references/metabolic_HMMs
assembly_list=~/Everglades/metadata/lists/2018_analysis_assembly_list.txt
ORFs=~/Everglades/dataEdited/assemblies/ORFs/

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $2 }' | while read HMM
do
  cd ~/Everglades/dataEdited/2018_analysis_assembly/
  geneName=$(awk -F ',' -v HMM="$HMM" '$2 == HMM { print $1 }' $metabolic_HMMs.csv)

  if [ ! -e metabolicProteins/identification/$geneName.afa ]; then
    echo "Searching for" $geneName
    mkdir metabolicProteins/identification/$geneName

    cat $assembly_list | while read assembly
    do
      if [ ! -e metabolicProteins/identification/$geneName/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout metabolicProteins/identification/$geneName/$assembly\_$geneName.out \
                  --cpu 4 \
                  --cut_tc \
                  $metabolic_HMMs/$HMM \
                  $ORFs/$assembly.faa \
                  > metabolicProteins/identification/$geneName/$assembly\_$geneName.txt
        lineCount=`wc -l < metabolicProteins/identification/$geneName/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $gene "hits in" $geneName
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  metabolicProteins/identification/$geneName/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  metabolicProteins/identification/$geneName/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM

    cd metabolicProteins/identification

    shopt -s nullglob
    for i in $geneName/*$geneName.faa; do FOUND=$i;break;done

    if [ ! -z $FOUND ]; then
      echo "Concatenating and aligning" $geneName
      cat $geneName/*$geneName.faa \
          > $geneName\_all.faa
      hmmalign -o $geneName.sto \
                  $metabolic_HMMs/$HMM \
                  $geneName\_all.faa
      $scripts/convert_stockhold_to_fasta.py $geneName.sto
      grep '>' $geneName\_all.faa | \
        sed 's/>//' \
        > $geneName\_all_list.txt
    else
      echo "No" $geneName "sequences found at all :("
    fi
    FOUND=""
  else
    echo "Already pulled out" $geneName "sequences"
  fi
done

conda deactivate

exit




#########################
# Extract depths of all scaffolds
#########################

screen -S EG_metabolic_genes_depth

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh

cd ~/Everglades/dataEdited/2018_analysis_assembly
mkdir metabolicProteins/depth

# Generate list of all scaffolds that have potential metabolic genes on them
grep '>' -h metabolicProteins/identification/*.afa | \
  sed 's/>//' | \
  cut -d"_" -f1,2 | \
  sort | \
  uniq \
  > metabolicProteins/depth/scaffold_all_list.txt


# Pull out all depths

cat ~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt | while read metagenome
do

  rm -f metabolicProteins/depth/$metagenome\_depth_raw.tsv
  conda activate bioinformatics
  PERL5LIB=""

  # Calculate depth over each residue in each scaffold
  cat metabolicProteins/depth/scaffold_all_list.txt | while read scaffold
  do
    assembly=$(echo $scaffold | awk -F '_' '{ print $1 }')
    echo "Calculating coverage of" $metagenome "over" $scaffold
    samtools depth -a -r $scaffold mapping/$metagenome\_to_$assembly.bam \
        >> metabolicProteins/depth/$metagenome\_depth_raw.tsv
  done
  conda deactivate

  # Average the depth of each residue over the entire
  echo "Aggregating" $scaffold "depth information for" $metagenome
  conda activate py_viz
  PYTHONPATH=""
  python ~/Everglades/code/generalUse/calculate_depth_contigs.py \
            metabolicProteins/depth/$metagenome\_depth_raw.tsv \
            150 \
            metabolicProteins/depth/$metagenome\_depth.tsv
  conda deactivate
  rm -f metabolicProteins/depth/$metagenome\_depth_raw.tsv
done



#########################
# Dereplicate all genes
#########################

screen -S EG_metabolic_genes_dereplication
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins
mkdir derep
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
metabolic_HMMs=~/Everglades/references/metabolic_HMMs
scripts=~/Everglades/code/generalUse

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $1 }' | while read gene
do

  if [ -e identification/$gene\_all.faa ]; then
    if [ ! -e derep/$gene\_clean.faa ]; then
      echo "Working on" $gene
      # Dereplicate sequences
      cdhit=~/programs/cdhit-master
      $cdhit/cd-hit -g 1 \
                    -i identification/$gene\_all.faa \
                    -o derep/$gene.faa \
                    -c 0.97 \
                    -n 5 \
                    -d 0
      $cdhit/clstr2txt.pl derep/$gene.faa.clstr > derep/$gene\_cluster.tsv

      # Generate list of dereplicated sequences
      grep '>' derep/$gene.faa | sed 's/>//' > derep/$gene\_derep_list.txt
      # Remove asterisk from sequences
      sed 's/*//g' derep/$gene.faa > derep/$gene\_clean.faa
      # Align dereplicated sequences to HMM
      HMM=$(awk -F ',' -v gene="$gene" '$1 == gene { print $2 }' $metabolic_HMMs.csv)
      hmmalign -o derep/$gene.sto \
                  $metabolic_HMMs/$HMM \
                  derep/$gene\_clean.faa
    else
      echo "Already dereplicated" $gene
    fi
  else
    echo "Nothing to dereplicate for" $gene
  fi
done

conda deactivate

# Convert sto format to fasta
cd ~/Everglades/dataEdited/2018_analysis_assembly/metabolicProteins/derep

metabolic_HMMs=~/Everglades/references/metabolic_HMMs
scripts=~/Everglades/code/generalUse

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $1 }' | while read gene
do
  if [ ! -e $gene.afa ]; then

    echo "Cleaning alignment for" $gene

    HMM=$(awk -F ',' -v gene="$gene" '$1 == gene { print $2 }' $metabolic_HMMs.csv)

    $scripts/convert_stockhold_to_fasta.py $gene.sto
    python $scripts/cleanFASTA.py $gene.afa
    mv -f $gene.afa_temp.fasta $gene.afa
    rm -f $gene.sto
  else
    echo $gene "alignment already converted"
  fi
done
