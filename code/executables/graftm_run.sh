#!/bin/sh

##############################
# code/executables/graftm_run.sh
# Benjamin D. Peterson

# This is a general script to run graftM,
# given a metagenome with paired reads and
# a graftM package.
##############################

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate graftm
PYTHONPATH=""
PERL5LIB=''

cd $workingDirectory

graftM graft --forward $metagenomeFolder/$metagenome\_R1.fastq.gz \
              --graftm_package $gpkg \
              --reverse $metagenomeFolder/$metagenome\_R2.fastq.gz \
              --input_sequence_type nucleotide \
              --log $metagenome.log \
              --threads 7 \
              --output_directory $metagenome\_output
