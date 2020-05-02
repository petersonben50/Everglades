### Protocol for assembly mapping to 2018 assemblies

This is the protocol I'm using for mapping reads to my assemblies from the 2018 dataset.

**Map reads and process output**

I used bowtie2 for my mapping.
First I had to build an index of the assembly.
Then I mapped the paired-end reads, the single reads, and the merged reads to the indices.
Once the mapping was done, I used samtools to convert the files to BAM files, then sorted and indexed the files, and finally removed the SAM and raw, unsorted BAM files.
