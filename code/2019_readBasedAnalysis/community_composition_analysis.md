### Workflow for community composition analysis

I wanted to capture the microbial diversity patterns in a more comprehensive manner than ones based on assemblies.
This workflow uses GraftM to identify and classify 16S reads in the metagenomes.


**Diversity analysis of metagenomes using Nonpareil**

I used Nonpareil 3 to estimate the coverage of the metagenomes and calculate the Nd, the Nonpareil diversity metric [@rodriguez-r2018].
I based my analyses on Johnston et al, 2019 in PNAS [@johnston2019].
I ran the program through a submission file, `nonpareil3_run.sub`.
Only the forward reads were used for analysis.
I used 100,000 reads for each metagenome.
The program was run on the kmer setting, with a kmer length of 32.
In the R script, I saved out the coverage of each metagenome into a table, plotted the Nonpareil curves all together, and plotted the Nonpareil diversity metrics for each metagenome.
