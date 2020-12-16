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
Generally, I set the number of clusters to somewhere between a half and a third of the predicted number of complete genomes.


**Add taxonomic information**

I also wanted to add information from GTDB to gain real-time insight into the taxonomy of the bins I was generating.
I did this in anvio, as outlined [here](http://merenlab.org/2019/10/08/anvio-scg-taxonomy/).
Make sure the GTDB databases are installed already.
Beyond that, just need to run `anvi-run-scg-taxonomy` on the database.


**Run CONCOCT binning**

We're going to focus in on using CONCOCT clustering to start the manual binning of the hgcA+ bins.
I'll use the Excel sheet I put together in the genome count estimation section to limit the number of clusters that we get.

*Troubleshooting*: I'm getting an error from anvio. I had not install concoct through conda, but just was using the default one from GLBRC.

*Side note*: Never use excel to make a csv that you're later going to use in linux.
Download Libre Office for this.


**Add MetaBat2 info**

Add the Metabat2 information as a collection.

**Search for hgcA+ bins**

Now, for the moment of truth: did any hgcA sequences make their way into any bins?
Let's find out.
First I summarized the bins and extracted a list of the original bin names.
I then searched for bins that had a list of hgcA in them, and saved the assembly ID and the bin name to a file (`original_hgcA_bin_list.txt`).
I downloaded this list to my local computer, and copied it into a new spreadsheet to use for my notes on the binning process: `dataEdited/2019_binning/binning_initial/anvioDB_processing/hgcA_binning_notes.xlsx`


**Manually bin hgcA+ bins**

Then I manually binned these!
Before I did, though, I copied the anvio database folder to a new folder (anvioDBs_modified) and modified the new one.
Then I went ahead and manually binned the hgcA+ bins.

Tried to manually bin them.
Overall, I'd say this process SUCKED.
There really aren't any definitive hgcA+ bins in here.
Nothing over 40% complete, and even to get there it was a little sketchy.
Might want to give up on the genome-resolved metagenomes dream for these metagenomes.
Maybe with the expanded number of metagenomes from our 2019 sampling, we'll be able to retrieve some, but I think we're stuck with assembly-based analyses for the 2019 work.
