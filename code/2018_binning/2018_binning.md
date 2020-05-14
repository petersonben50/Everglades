### Protocol for initial binning workflow on 2018 EG metagenomes

**Prepare scaffolds and mapping files**

*Filter out short scaffolds*

First, I used anvio to filter out any scaffolds that were shorter than 2000bp.

*Map reads to filtered scaffolds*

I then mapped all the reads from the metagenomes back to these subsetted scaffolds using bowtie2 (v2.2.2).
First, I built the indices for each assembly, then mapped all the forward, reverse, single, and merged reads to each index.
I then used samtools to convert sam files to bam files, which I then sorted and indexed.

**Run automatic binning algorithms**

I then wanted to use an automatic binning method to generate bins, both for their own use and for help in guiding the manual anvio binning.
I started with Metabat2, and will only use this one for now.
First, I used the JGI script to summarize the bam files for all metagenomic mapping to each assembly into a depth file.
I then binned the scaffolds using default settings.
Finally, I renamed the bins to include their assembly of origin and metabat.

**Generate anvio databases**

For the hgcA+ binning, I wanted to do it manually in anvio.

*Generate contig databases*

The first thing I needed to do was generate the contig databases, which I did using anvio v6.2 (I'll use that version for the rest of the binning as well).
I generated the contig database with default settings.

*Populate contigs DB with HMMs*

I then populated the DB with different HMMs.
I used the default single-copy gene HMM set in anvio that I'll use for completeness and redundancy estimates.
I also used the HMM folder that I built for HgcA to populate the database with hgcA hits.

*Add kaiju annotation information*

I also wanted to add in taxonomic information for each of the scaffolds using Kaiju, basically following [this workflow](http://merenlab.org/2016/06/18/importing-taxonomy/#kaiju).
Briefly, I used anvio to retrieve nucleic acid sequences for gene calls.
I then classified each scaffold with kaiju using the NCBI nr database downloaded on March 31st, 2020.
I then added this information back into anvio.

**Generate read profiles**

Then, I profiled the reads mapping to each assembly and added that to the anvio database.
When merging the reads, I skipped the hierarchical clustering and the CONCOCT binning.
I'll perform binning with a condensed number of clusters later on.

**Estimate number of genomes**

Before performing our clustering, I wanted to get an estimate of the number of genomes that we had in each assembly.
I used the anvi-display-contigs-stats program to do this, then manually saved the estimated number of genomes to `estimated_number_of_genomes.csv` (in the second column).
In the third column, I added the number of clusters I'd like CONCOCT to generate for me.
