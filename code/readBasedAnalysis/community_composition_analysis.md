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


**Mash analysis of beta diversity**

I calculated the mash distances in the `assembly_processing` scripts.
Then I ordinated these distances (using PCoA) in `mash_PCoA_seds.R`.
Finally, I wanted to statistically test the distances between these communities, preferably using a PERMANOVA.
They are statistically different, obviously, at the different sites.
Could probably test this against the sulfide as well.
I wonder though, because we didn't measure it in each of the cores and instead have a single sulfide variable for each site, is that really a valid computation, since we're artificially adding confidence to the link between the microbial community and the sulfide levels.
Hmm, with this dataset, since the metagenomes are technical replicates, we'd probably need to average the data together, maybe?
I don't know, need to ask this of a real ecologist.
 
