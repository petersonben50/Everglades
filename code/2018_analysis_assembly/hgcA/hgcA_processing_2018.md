### My workflow for analyzing the hgcA sequences in the Hell's Canyon assemblies

#### hgcA identification


**Identify potential hgcA seqs**

I searched for hgcA sequences in the assemblies using my hgcA HMM, through hmmsearch, using the trusted cutoff built into the HMM.
I extracted the amino acid sequences that hit the HMM using a python script and concatenated them.
Final count of 50 hgcA sequences, but many of these are going to be the same ones.

**Align and quality filter results**

I aligned all the putative hgcA sequences to the hgcA HMM, converted it to fasta format, and downloaded it to my computer for manual inspection.
I manually inspected the sequences in Geneious, looking first for the cap helix domain, then for the transmembrane domains (at least 2).
The MegaHit coassembly churned out 8 partially assembled hgcA sequences.
Wonder if it's due to the increased strain variation leading to fragmentation?
The Pw04 MegaHit assembly returned 3 fragments, one from Pw03 MegaHit, and 1 each from Pw01 and Pw02 in the metaSPADes assemblies.
There's one sequences (Pw04Meta18_000000032264_1) that only has one predicted TM domain, and doesn't align well with the others, but isn't cut off at the C-terminus.
I decided to hang on to that one for now.
There is a scaffold (Pw01Meta18_000000017085) with two hgcA hits against it, at the first and second open reading frame.
It appears that the ORF was split, likely a Prodigal artifact.
I'm just going to cut these two out rather than trying to combine them.
During the depth analysis (below) it turned out that the scaffold wasn't at very high coverage anyways.

I saved a file (`hgcA_good.afa`) with all the keeper sequences, and only kept those with a full cap helix domain and several transmembrane domains.



**Extract gene neighborhood data for hgcA+ scaffolds**

I also wanted to inspect the position and gene neighborhood of the hgcA sequences on their scaffolds.
To do this, I first pulled out the scaffolds and the corresponding GFF entries for each hgcA+ scaffold.
I then isolated the gene neighborhoods using a python script I modified from Tyler Barnes.
I pulled out sequences 5000bp from before and after the hgcA sequence, as well as the corresponding GFF entries.
I did this for both the raw list of hgcA sequences as well as the trimmed list.



**Search for downstream hgcB sequences**

First I pulled out the genes downstream from hgcA using a custom python script (`retrieve_downstream_gene_name.py`).
I searched the corresponding downstream genes using a hgcB HMM I built for the 5M project using the hgcB recovered from my assemblies.
After using a few different cutoffs and checking the results, I ended up using 30 as a cutoff score.
I aligned the hits against the HMM, and they all had the CM/IECGAC region.
I saved out a list of the hgcB names.


**Extract depths of hgcA**

Next I wanted to extract the coverage information for all these sequences.
I used the depth function in samtools to calculate the coverage over each nucleotide residue of the hgcA+ scaffolds.
I then used the `calculate_depth_contigs.py` script to calculate the average coverage of each nucleotide over the sequence, excluding 150 bp on both ends of the scaffold, to eliminate any bias against shorter contigs (in case reads aren't mapping to the ends of the scaffold as well).
I did this for each sequence initially identified as hgcA, so that I could see if any of the incomplete sequences are relatively abundant.

I downloaded the depth files to my local computer, then combined them into a single file and normalized them with the R script `clean_hgcA_abundance.R`.
This saves out a file with the aggregated and normalized depth information here: `dataEdited/2018_analysis_assembly/hgcA/depth/hgcA_coverage.csv`.


**Clustering of sequences**

Next we'll want to dereplicate the hgcA sequences across our assemblies.
We'll use a cutoff of 97% identity for this process.
In addition to getting the fna file, we'll generate the text file with the cluster information, then download that to my computer.
I also did this using a cut-off of 80% identity to look at patterns of sequence similarity.

**Dereplicate sequences**

Next I wanted to pull all the data I had on hgcA sequences into R to get a good look at them.
I used an Rmd file to keep track of my notes on this: `hgcA_dereplication.Rmd`.
Bottom line: we'll go with the hgcA representatives that CD-HIT returned.


**Make phylogenetic tree**

*Use McDaniel et al 2020 references*

First I wanted to flesh out our hgcA tree using some reference sequences.
I downloaded the hgcA.faa file from Elizabeth's Hg paper, and converted the file to a blast database.
I then blasted our sequences against that database, keeping the top 5 hits.
I concatenated our sequences with the references that we found.
With this set, I generated an alignment using MUSCLE (v3.8.31).
To check the dataset, I generated a rough tree using FastTree (v2.1.10) to check on the quality of the alignment.
I downloaded this to Geneious and checked it out.
Looking a little thin with just the references from Elizabeth's paper, so I decided to check it out in R to see what it needs.
Code is here: `code/2018_analysis_assembly/hgcA/clean_hgcA_tree.R`.

The branches were a little long for the sequences near Deltaproteobacteria, and the one clustering with Elusimicrobia.
So, I decided to flesh it out with the nr database

*Use nr database to find references and confirmed methylators*

To flesh out the tree, I blasted my sequences against the non-redundant NCBI database (downloaded 2020-05-07), keeping the top 5 results.
I dereplicated these sequences against the sequences from McDaniel et al 2020, pulled out the unique sequences.
Finally, I concatenated the McDaniel 2020, the nr sequences, and the 30 confirmed methylators into a single file along with the assembly sequences.
I aligned them using MUSCLE, then generated a tree using FastTree.
I downloaded the tree file to my local computer, where I visualized it in the `clean_hgcA_tree.R` code file.

It was overpopulated at this point.
I manually inspected the tree and identified the accession numbers of the files I wanted removed, which I saved here: `dataEdited/2018_analysis_assembly/hgcA/phylogeny/seqs_to_remove_tree_2.txt`.


*Generate final tree*

I then generated the final list of sequences to use in my tree file.
I first removed all the sequences that I didn't want from the second iteration (saved in `seqs_to_remove_tree_2.txt`), using a python script I wrote: `code/generalUse/remove_fasta_seqs_using_list_of_headers.py`.
I then also added in the hgcA sequences from the SYN bins that I found in my Mendota metagenomes, since the Syntrophobacterales branches seemed to be a little thin.
Then I generated an alignment using MUSCLE and made a tree using FastTree, just to check it out.
Checked it out in the `clean_hgcA_tree.R` file.
This looks pretty damn good actually, now just need to clean up the alignment before feeding it to RAxML.

To clean the alignment, I masked it with with BMGE1.1,
I did this on the [Galaxy portal](https://galaxy.pasteur.fr/?tool_id=toolshed.pasteur.fr%2Frepos%2Fdcorreia%2Fbmge%2Fbmge%2F1.12&version=1.12&__identifer=5anq311b6rv).
I used the BLOSUM30 substitution matrix, with everything else on default (sliding window of 3, a maximum entropy threshold and a gap rate cut-off of 0.5, and a minimum block size of 5).
I then manually checked the alignment in Geneious.
Several sequences were truncated on either end, so I cut 48 residues off the front end and 77 off the back of the alignment.
The final exported sequence is `hgcA_for_phylogeny_final_trimmed_cut.afa`.


*Generate tree in RAxML*

I then generated an ML tree using RAxML.
We used RAxML (v8.2.11) for this.
The seed value was 283976
Rapid bootstrapping was used (with the "-N autoMRE" flag, and a seed number of 2381).
Automatic detection of convergence stopped bootstrapping after xxx replicates.
The phylogeny was generated under a gamma distribution and the best substitution matrix was automatically determined to be LG (-m PROTGAMMAAUTO).

*Clean tree*

The tree was then downloaded and loaded into R (`clean_hgcA_tree.R`).
It was mid-point rooted using the phangorn package and visualized using ggtree.
I replaced the accession numbers of the references with their names, and saved out a cleaned hgcA tree (`dataEdited/2018_analysis_assembly/hgcA/phylogeny/hgcA_clean_tree.rds`).
I'll also use this script to save out a color vector for the tree.
I used these two RDS files, along with the depth information for the hgcA+ scaffolds, to generate an hgcA tree paired with the overall abundance information.


**hgcA phylogenetic cluster assignment**

I then assigned a taxonomic group to each of the hgcA sequences, based on the phylogenetic tree that I had built.
That information was saved here: `hgcA_phylogenetic_clusters.xlsx`.
I started out by classifying them simply as methanogen, syntroph, Geobacter, or Elusimicrobia.
At some point I might want to split up the methanogens into smaller categories.




**hgcA scaffold analysis**

I set up an analysis script in `hgcA_scaffolds.Rmd` to look at the scaffold coverage and length.
Check out my notes there.
