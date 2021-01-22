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
Also added a script to run an analysis of variance on the diversity metrics.
Nothing here.


**16S analysis of metagenomes using GraftM**

*Run GraftM*

I ran GraftM through submission files on GLBRC, using GraftM version 0.13.1.
I used the Silva v132 database 7.71.silva_v132_alpha1.gpkg as a reference database for this.
I ran `graftM graft` using the default setting.

*Clean up GraftM outputs*

I copied out the `combined_count_table.txt` file.
Download all these files to my local computer: `dataEdited/2019_readBasedAnalysis/16S/graftM`
Clean them using R: `code/2019_readBasedAnalysis/cleaning_graftm_data.R`.
Finally, I generated a plot of the phyla present here: `code/2019_readBasedAnalysis/graftm_phyla_plots.R`.
