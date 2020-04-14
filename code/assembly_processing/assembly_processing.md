### Assembly Processing for 2019 intensive

This document outlines the protocol I followed for processing the metagenomes and assemblies from the Everglades.
It will also include my notes about the process, and directions to other files as needed.

**Data retrieval and inspection**

First thing I did was retrieve the data.
The first round of sequencing was through Michigan State's sequencing center.
I was able to retrieve it using `wget`.
All these metagenomes are in `~/Everglades/dataRaw/metagenomes`.

I next made a list with all the metagenomes of interest: `metadata/lists/metagenome_list.csv`.
This list will be used to check and trim the metagenomes.


**Trimming metagenomes**

Next I quality trimmed the metagenome.
To do this, I used the fastp program (version 0.20.0), which trims the metagenome to our liking and gives us the stats about the metagenome both before and after trimming.
I download the program using conda, into our bioinformatics virtual environment.

I used this program to save out both the paired and unpaired reads after trimming.
I also stored the failed reads, just so we can check to see what's driving the loss of reads if need be.

I cut the tail end of the reads based on the read quality.
The front end quality of the reads are usually pretty good so I didn't bother with trimming the front or using the `cut_right option`.
I used a 10bp window with a quality score cutoff of 20.
Since I was trimming these down, I wanted to cut out any short reads that get overtrimmed.
Length filtering is set to 15 by default, but I changed it to 100.
I left the default quality filtering in place as well, which is a quality limit of 15, and a maximum of 40% of bases in a read below that threshold.
I also merged the reads.
For some reason, once I use the merge function, the single reads are filled in at the end, so don't get worried if you're checking the size of the files as the program runs and the single reads file is empty.
I then went through and took a look at the output of fastp, the files of which are stored here: `/Users/benjaminpeterson/Documents/research/Everglades/dataRaw/metagenomes/2019_intensive/fastp`.

-




Unfortunately fastp doesn't read out how many single reads vs. paired end reads there are, so I ended up running fastQC on the data as well just to get those counts and get some other information on the read quality.
I'll also save the read counts here: `~/Everglades/dataEdited/metagenomes/2019_intensive_ancillary_info/metagenome_read_count.tsv`, then download it to `dataEdited/metagenomes/reports/metagenome_read_count.tsv`.






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

I didn't do the surface water samples from 2018 for now, since in our initial analysis they weren't of much interest.

Also, before starting assemblies, I needed to upload the two `assembly_by_group` python scripts for MegaHit and metaSPADes, to `code/assembly_processing`.

*Assemblies by metaSPades*

Once the metagenomes were trimmed and ready to go, I assembled them individually using metaSPADes (v3.13.2).
For now, I focused on getting individual assemblies with metaSPADes, since the coassembly was running out of memory.
I included the paired reads, the merged reads, and the single reads.
For kmers, I included 21, 33, 55, 77, 99, and 127.

*Assemblies by MegaHit*

I also assembled the metagenomes using MegaHit.
For this, I used MegaHit version 1.2.9.
I assembled each of the 2018 porewater samples individually and together as a coassembly.
I used kmers from 21 to 121 with a step of 10.


**Clean up metagenome assemblies**

For each of my assemblies, I used anvio to keep only scaffolds that were longer than 1000bp and renamed the scaffolds to include the assembly ID in the scaffold name.
This was done in separate for loops for MegaHit and metaSPADes.
