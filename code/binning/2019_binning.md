### Protocol for initial binning workflow on 2019 EG metagenomes

**Prepare scaffolds and mapping files**

*Filter out short scaffolds*

First, I used anvio to filter out any scaffolds that were shorter than 2000bp.
I ran this in a for loop, since it doesn't take too long.

*Map reads to filtered scaffolds*

I then mapped all the reads from the metagenomes back to these subsetted scaffolds using bowtie2 (v2.2.2).
I ran this as a submission job, since it takes so long to map things.
Set it up to run the mapping to each assembly in a different job, which allowed me to set up the indices within each run as well.
In the executable files, I first built the indices for each assembly, then mapped all the forward, reverse, single, and merged reads to each index.
I then used samtools to convert sam files to bam files, which I then sorted and indexed.


**Generate anvio databases**

For the hgcA+ binning, I wanted to do it manually in anvio.
I started generating these databased in parallel to doing the mapping.
These script are all run as a for loop, since they don't really take that long.

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



**Run automatic binning algorithms**

I then wanted to use an automatic binning method to generate bins, both for their own use and for help in guiding the manual anvio binning.
I started with Metabat2, and will only use this one for now.
First, I used the JGI script to summarize the bam files for all metagenomic mapping to each assembly into a depth file.
I then binned the scaffolds using default settings.
The bins are named to include their assembly of origin and metabat.
Finally, I removed the periods from the names and replaced them with underscores.




**Generate read profiles**

Then, I profiled the reads mapping to each assembly and added that to the anvio database.
This is done in a submission script, where the executable is generalized and the submission file can be customized to different years or projects.
When merging the reads, I skipped the hierarchical clustering and the CONCOCT binning.
I'll perform binning with a condensed number of clusters later on.



**Estimate number of genomes**

Before performing our clustering, I wanted to get an estimate of the number of genomes that we had in each assembly.
I used the [anvi-display-contigs-stats](http://merenlab.org/2016/06/22/anvio-tutorial-v2/#anvi-display-contigs-stats) program to do this, then manually saved the estimated number of genomes to `estimated_number_of_genomes.csv` (in the second column).
In the third column, I added the number of clusters I'd like CONCOCT to generate for me.
There's an Excel sheet too, with more details, like how many bacteria vs. archaea were predicted: `estimated_number_of_genomes.xlsx`.
I divided the estimated number of clusters by 3 and used that as the maximum number of CONCOCT clusters requested.


**Add taxonomic information**

I also wanted to add information from GTDB to gain real-time insight into the taxonomy of the bins I was generating.
I did this in anvio, as outlined [here](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/).
Make sure the GTDB databases are installed already.
Beyond that, just need to run `anvi-run-scg-taxonomy` on the database.


**Run CONCOCT binning**

We're going to focus in on using CONCOCT clustering to start the manual binning of the hgcA+ bins.
I'll use the Excel sheet I put together in the genome count estimation section to limit the number of clusters that we get.
This is run as a for loop.


**Add MetaBat2 info**

Add the Metabat2 information as a collection.

**Search for hgcA+ bins**

Now, for the (first) moment of truth: did any hgcA sequences make their way into any bins?
Let's find out.
First I summarized the bins and extracted a list of the original bin names.
I then searched for bins that had a list of hgcA in them, and saved the assembly ID and the bin name to a file (`original_hgcA_bin_list.txt`).
Wow. We have a ton here!
I downloaded this list to my local computer, and copied it into a new spreadsheet to use for my notes on the binning process: `dataEdited/2019_binning/binning_initial/anvioDB_processing/hgcA_binning_notes.xlsx`


**Manually bin hgcA+ bins**

Then I manually binned these!
(Second moment of truth).
Before I did, though, I copied the anvio database folder to a new folder (anvioDBs_modified) and modified the new one.
Then I went ahead and manually binned the hgcA+ bins.
Notes can be found here:
`/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/2019_binning/binning_initial/anvioDB_processing/hgcA_binning_notes.xlsx`.


**Summarize/export curated bins**

*Rename and summarize bins*

Fist, I renamed the bins to include the assembly ID.
These new bins were saved to the refined_and_renamed collection.
I then summarized the bins to get the output we needed.

*Pull out DNA files from hgcA+ bins*

I only pulled out the hgcA+ bins from our binning.
There are a lot of junky ones in here too, but just going to pull everything out for now.



**Check quality of bins**

I checked the quality of these bins using two metrics, the completeness/redundancy metrics from anvi'o as well as those from CheckM.

*Completeness/redundancy estimates from anvio*

I concatenated the summaries into a single folder, and then into a single file.
I then made a folder with the stats for just the hgcA+ bins.
Finally, I made a separate summary file and a list for just the hgcA+ bins with completion scores > 50 and redundancy < 10.

*Completeness/redundancy estimates from CheckM*

I also ran checkM on the bins that I generated. Used the `lineage_wf` workflow.
Then I went for the checkM qa script to save out the output file.


**Check out taxonomy of bins with GTDB**

I used the GTDBTK `classify_wf` program to classify the good bins, and summarized the output file.
I then downloaded the `taxonomy_summary.txt` file and manually generated an excel sheet with a taxonomic classification for each bin: `taxonomy_summary.xlsx`.
We'll use this to attach a name to binned hgcA sequences in the assembly-based tree.

**Compare bins by differential coverage**

To make sure that the high matching sets I generated down the line were appropriate, I extracted the coverage of all the bins from anvio for use in grouping bins by differential coverage.
I saved out a depth file with only the coverage of the good hgcA+ bins I was including.

**Get ORFs for bins**

I used Prodigal to predict the open reading frames for each of the bins.
I don't like extracting the ORFs from anvi'o, probably because those are run on the metagenomic mode, I think.
I run Prodigal on the single genome mode for these purposes, and save out the faa, fna, and GFF outputs.

**Run ANI comparisons on good bins**

I ran Sarah Steven's ANI [DAGman calculator](https://github.com/sstevens2/ani_compare_dag) on the good hgcA+ bins to generate the high-matching sets.

**Dereplicate bins**

I then aggregated the bin information in R, here: `code/2019_binning/binning_stats.R`.
I saved out a CSV file with all the bin information and manually dereplicated the hgcA+ bins.
Once I have this final list, I saved out the list here: `dataEdited/2019_binning/binning_initial/binsFinal_list.txt`.
I uploaded this file to GLBRC: `~/Everglades/dataEdited/2019_binning/binning_initial/binsFinal_list.txt`.
Just to make everything easier, I copied the scaffolds and the ORFs for these bins into a new `binsFinal` folder.
I also built a scaffold-to-bins file and a genes-to-bins file.

**Confirm *hgcA* in bins**

First I generated a key to relate the hgcA sequences from the bins to the seqs from the assembly-based methods.

*hgcA search*

I then manually searched the bins as well, just to confirm the quality of the hgcA sequences.
I also read the key into the R script I have to include the bin names in the assembly-based hgcA tree.
I did some manual analysis on these, saved out here: `dataEdited/2019_binning/hgcA/notes_on_bin_hgcA_seqs.xlsx`
In the HMM analysis, Sed991Mega19_000001090714_3 is showing a very low score (137), barely over the threshold.
This is one of the fused sequences
I also pulled out the sequences and aligned them to the HMM.
Sed993Mega19_000002407322_10 and Sed993Mega19_000002321000_3 both have the cysteine residue replaced with a serine group.
Interestingly, these are in the two Spirochaetes bins.
These are two of the four *hgcA* paralogs that I identified in the assemblies.

The hgcA sequence from Sed992Meta19_Bin_00002 (Sed992Meta19_000000003684_8) is truncated, and I didn't include it in the original assembly.
I'm going to include the bin in our analysis, for completeness sake, but this one won't show up in the hgcA tree.


**Run metabolic HMMs on bins**

I used my batch HMMs scripts to run a whole set of metabolic HMMs on the ORFs from the bins.
I also ran Shaomei's MHC python script on the ORFs.
