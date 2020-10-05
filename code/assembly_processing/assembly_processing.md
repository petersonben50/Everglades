### Assembly Processing for 2019 intensive

This document outlines the protocol I followed for processing the metagenomes and assemblies from the Everglades.
It will also include my notes about the process, and directions to other files as needed.

**Data retrieval and inspection**

First thing I did was retrieve the data.

*2018 data*:
The first round of sequencing was through Michigan State's sequencing center.
I was able to retrieve it using `wget`.
All these metagenomes are in `~/Everglades/dataRaw/metagenomes`.

I next made a list with all the metagenomes of interest: `metadata/lists/metagenome_list.csv`.
This list I made manually.
This list will be used to check and trim the metagenomes.

*2019 data*:
This sequencing was all done at QB3.
I retrieved it using lftp.
I then added these metagenome IDs to `metadata/lists/metagenome_list.csv`.
For the KMBP005 samples, they added an underscore for some reason.
Let's remove that.


**Trimming metagenomes**

Next I quality trimmed the metagenomes.
To do this, I used the fastp program (version 0.20.0), which trims the metagenome to our liking and gives us the stats about the metagenome both before and after trimming.
I download the program using conda, into our bioinformatics virtual environment.

I used this program to save out the paired reads, unpaired reads, and merged reads after trimming.
I also stored the failed reads, just so we can check to see what's driving the loss of reads if need be.

I cut the tail end of the reads based on the read quality.
The front end quality of the reads are usually pretty good so I didn't bother with trimming the front or using the `cut_tail` option.
I used a 10bp window with a quality score cutoff of 20.
Since I was trimming these down, I wanted to cut out any short reads that get overtrimmed.
Length filtering is set to 15 by default, but I changed it to 100.
I left the default quality filtering in place as well, which is a quality limit of 15, and a maximum of 40% of bases in a read below that threshold.
I also merged the reads.
For some reason, once I use the merge function, the single reads are filled in at the end, so don't get worried if you're checking the size of the files as the program runs and the single reads file is empty.
I then went through and took a look at the output of fastp, the files of which are stored here: `/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/metagenomes/reports/fastp`.

Unfortunately fastp doesn't read out how many single reads vs. paired end reads there are, so I also ran fastQC on the data as well just to get those counts and get some other information on the read quality.
Data is stored locally here: `/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/metagenomes/reports/fastqc`.

**Check size of metagenomes**

I also got a read count by counting headers in the trimmed fastq files, and saved it to `/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/metagenomes/reports/metagenome_read_count.tsv`.
I did this for the pre- and post-trimmed metagenomes.

Finally, I wanted to calculate the total number of nucleotides from each metagenome, which I could use to normalize the coverage of all our sequences.
For this, I used the [readfq program](https://github.com/billzt/readfq).
For the pre-trimming reads, I counted the coverage in both the forward and reverse reads.
For the post-trimming reads, I counted the number of nucleotides in the forward, reverse, single, and merged read files, all of which will be used for mapping.
I downloaded all these files to my local computer (`/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/metagenomes/reports/metagenome_coverage.tsv`).

In the R script `code/assembly_processing/mapping_normalization.R`, I generated an R object containing a vector that has normalization values for each metagenome.
I did this by normalizing by coverage: since I'm generating my abundance values by coverage rather than read mapping, I should be normalizing by the nucleotide residues rather than the number of reads.
So, I normalized the coverage to a 2x150bp metagenome with one million reads, for a normalized metagenome coverage of 3e+10.
I divided the total coverage of each of our metagenomes (from R1, R2, single, and merged) by this normalized coverage.

**Notes on metagenome quality**

First, I wanted to get a general sense of how many reads and how much coverage was trimmed out by fastp.
I checked that out in this R script: `code/assembly_processing/extent_of_trimming.R`.
The trimming cut out between 3 and 6% of the total coverage of the difference metagenomes.
None of them jump out as being anomalously high or low.


Then I went about checking the status of each metagenome, mostly by looking at the .
For the 2018 samples, I only looked at the porewater samples, not the surface water samples, since I'm not planning to use those for anything yet.
I looked at the fastQC, the fastp, the read counts, and the coverage data for this.

- ENP18_001_002_003: Looks pretty good. Single reads are of much lower quality than the forward and reverse ones. Probably to be expected. Fastp documents are much less intuitive than the fastQC ones.

- ENP18_024_025: Pretty similar to ENP18_001_002_003.

- ENP18_030_032: Nothing jumps out. Similar to previous ones.

- ENP18_048_049_50: Nothing really jumps out.

- ENP18_061: Seems to have slightly higher deviation from the expected GC distribution.

General takeaways:
The fastp output isn't super helpful, I like the fastQC data much more.
Single reads are generally of much lower quality than R1 and R2.
The GC distribution of the reads is not in line with the theoretical GC content.
There's a small enrichment of reads 138/139bp long, but the majority of them are definitely 150bp.
Overall, things look pretty good!
Happy to go ahead with the assemblies from here.


**Metagenome assembly**

I ran individual assemblies as well as a co-assembly for these metagenomes.
All using both metaSPADes and MegaHit, and will evaluate them afterwards to see which ones were the best.
To set this up, I made a spreadsheet for each assembly program to use (`metadata/lists/assembly_groups_metaspades.csv` and `metadata/lists/assembly_groups_megahit.csv`).
This spreadsheet has two columns, "assemblyID" and "metagenomeID".
Each metagenome that is to be a part of an assembly gets it's own row:
So, if coassembly_2 is going to use ENP18_1 and ENP18_2, this would is what the sheet would have:

assemblyID,metagenomeID
coassembly_2,ENP18_1
coassembly_2,ENP18_2

Also, before starting assemblies, I needed to upload the two `assembly_by_group` python scripts for MegaHit and metaSPADes, to `code/assembly_processing`.

**2018**: I didn't do the surface water samples from 2018 for now, since in our initial analysis they weren't of much interest.

**2019**: Just wanted to jot some notes down here on how I'm going to process the 2019 data.
I'm thinking to just assemble the sediment samples first, since I think those will be of greater interest.
We did do two metagenomes at each site, and will probably will be best to co-assemble those two metagenomes for each site.
Won't do any co-assemblies with all the metagenomes at this point, since the samples are so separated.
I'm also not planning to individually assembly each metagenome at this point.


*Assemblies by metaSPades*

**2018**: Once the metagenomes were trimmed and ready to go, I assembled them individually using metaSPADes (v3.13.2).
For now, I focused on getting individual assemblies with metaSPADes, since the coassembly was running out of memory.
I included the paired reads, the merged reads, and the single reads.
For kmers, I included 21, 33, 55, 77, 99, and 127.

**2019**: Focused on coassembling the two metagenomes from the same core, which essentially are technical replicates.
Same thing as 2018: paired, unpaired, and merged reads; 21, 33, 55, 77, 99, and 127 for kmers.

*Assemblies by MegaHit*

**2018**: I also assembled the metagenomes using MegaHit.
For this, I used MegaHit version 1.2.9.
I assembled each of the 2018 porewater samples individually and together as a coassembly.
I used kmers from 21 to 121 with a step of 10.

**2019**: Again, focused on assembling the two metagenomes from the same core.
Used the same method as in 2018, kmers 21 to 121 by 10, with the paired ends and the single reads (no merged reads here).


**Clean up metagenome assemblies, get stats**

For each of my assemblies, I used anvi'o to keep only scaffolds that were longer than 1000bp and renamed the scaffolds to include the assembly ID in the scaffold name.
This was done in separate for loops for MegaHit and metaSPADes.

I then got some stats on the size and quality of the assemblies.
For this, I used the `abyss-fac.pl` script included with the Abyss assembler.
I aggregated the numbers into `all_assemblies_stats.txt` and downloaded it to my local computer here: `dataEdited/assemblies/reports`.
Let's take a look.

**2018**: Generally, the MegaHit assemblies had a higher N50 than the metaSPADes ones.
However, the metaSPADes assemblies had a higher overall length.
For some reason, despite using anvi'o to trim out any scaffold under 1000bp, there are a few scaffolds in the metaSPADes assemblies that are under 1000bp, down to about 950bp.
I think I'm going to run my hgcA and other assembly-based analyses on all these assemblies and can then simply dereplicate across them.
I'll also do this for the binning process.

**Predict open reading frames**

I used Prodigal on the `meta` setting to predict the open reading frames.
I saved out the nucleic acid sequences, the amino acid sequences, and a GFF file, then cleaned up the fasta files using `cleanFASTA.py`.
Finally, I counted the number of open reading frames in each assembly, and downloaded the counts here: `dataEdited/assemblies/reports/ORF_counts.tsv`.
Generally, the metaSPADes assemblies have many more open reading frames than the MegaHit assemblies.
