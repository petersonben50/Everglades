### Protocol for assembly mapping to 2019 assemblies

This is the protocol I'm using for mapping reads to my assemblies from the 2019 dataset.

**Map reads and process output**

I used bowtie2 for my mapping.
First I had to build an index of the assembly.
Then I mapped the paired-end reads, the single reads, and the merged reads to the indices.
Once the mapping was done, I used samtools to convert the files to BAM files, then sorted and indexed the files, and finally removed the SAM and raw, unsorted BAM files.


**Calculate fraction of reads mapped to each metagenome**

I also wanted to know what fraction of each metagenome was mapping to each assembly.
For this, I used the `flagstat` function in `samtools` and extracted the mapped percentage out of the fourth row.
These were all saved in a long tsv file: `mapping_percentages.tsv`
I downloaded these to my computer, here: `/Users/benjaminpeterson/Documents/research/Everglades/dataEdited/2019_analysis_assembly/reports`.

**Calculate overall fraction of reads in a metagenome mapping to any assembly**

I also wanted to know the fraction of reads in a metagenome that mapped to any assembly, to get a sense for how much more coverage we're getting using the multi-metagenome approach.
I had to rig this analysis up a little.
First we'll read out the names of each read that was mapped to any assembly.
I used `samtools view` to do this, with the flag `-F 4` to remove any reads with a 4 bit in the flag field, which indicates an unmapped read.
However, the forward and reverse reads have the same read ID, so I needed to distinguish them.
I had to count them all individually.
I used an online database to determine the appropriate flags: https://broadinstitute.github.io/picard/explain-flags.html.
First, I did single reads (which will include merged reads), which requires excluding the flags 1 (paired read) and 4 (read unmapped).
So, flag is `-F 5`.
Next I pulled out the mapped reads that were the first in a pair.
To do this, I set that the read needed to be paired and be the first in the pair (`-f 65`) but that it could not be unmapped (`-F 4`).
For the reverse reads, required to be paired and be the second in the pair (`-f 129`) but that it could not be unmapped (`-F 4`).
And, just as a sanity check, let's pull out the unmapped reads into a separate file (`-f 4`).
The sum length of these four files should be the sum of the reads in the metagenome.

Then, within each type of file, I concatenated the unique mapped read names for each metagenome into a separate file.
I then counted all the lines in each type of file for each metagenome and saved those counts in the mapping reports folder: `uniq_mapped_reads_MG.tsv`.
I downloaded these to my computer, here: `dataEdited/2019_analysis_assembly/reports`.

I then moved into R to analyze these: `code/2019_analysis_assembly/mapping_stats.R`.
I read in the unique mapped reads count and the total reads counts from each metagenome.
We won't worry about the mapping percentages for now.
Not doing assemblies on the porewater is obvious, as only 2-13% of the porewater metagenomes map to any of the assemblies.
Wonder if it's worth redoing this analysis with the porewater assemblies? (No, not at this point, is the answer).
For the sediment metagenomes, we see overall mapping percentages from 19-44%.
2A-N and 3A-N have the highest percentage of mapped reads, while 2A-A and 3A-F seem to be the worst.
The good news is, is that this doesn't correlate directly to the higher levels of *hgcA* and metabolic proteins we see in some sites.
Not that I'd really expect that, but good to know anyways.
