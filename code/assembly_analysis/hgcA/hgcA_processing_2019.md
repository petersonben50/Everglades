### My workflow for analyzing the hgcA sequences in the Everglades assemblies

I went ahead with using both metaSPADes and MegaHit assemblies for each of the six sights.
I'd like to do a comparison in here to see if this is really necessary, or if MegaHit or metaSPADes pulls out more genes, but not going to prioritize that at this point.

#### hgcA identification


**Identify potential hgcA seqs**

I searched for hgcA sequences in the assemblies using my hgcA HMM, through hmmsearch, using the trusted cutoff built into the HMM.
I extracted the amino acid sequences that hit the HMM using a python script and concatenated them.
Final count of 50 hgcA sequences, but many of these are going to be the same ones.

**Align and quality filter results**

I aligned all the putative hgcA sequences (`hgcA_raw.afa`) to the hgcA HMM, converted it to fasta format, and downloaded it to my computer for manual inspection.
I manually inspected the sequences in Geneious, looking first for the cap helix domain, then for the transmembrane domains (at least 2).
- 991 Mega: 3 seqs original, two seem to be the weird divergent ones. Will keep all though.
- 991 Meta: 7 originally. Two have ends truncated, one has front truncated. Two are those divergent seqs.
- 992 Mega: 13 sequences. Lot of truncated ones here. 4 are truncated badly from end, one truncated from the front. 3 more have minor truncation at end, but enough to kick them out. Left with 5.
- 992 Meta: 11 sequences. Better quality. Only 4 are truncated, all on C-terminus end
- 993 Mega: 58 seqs here. Many incomplete ones. Cut out 13 for not having TM regions, they were incomplete. Should look to see how many of these were predicted to be truncated by Prodigal. 6 more were completely missing cap helix domain. Three more of these sequences (Sed993Mega19_000002407322_10, Sed993Mega19_000002321000_3, Sed993Mega19_000001097352_5) were divergent at that cap helix domain, but were pretty homologous otherwise, and had the predicted TM domains. So, we'll include them for now, but should keep an eye on them. I want to include them in the tree. 39 total sequences
- 993 Meta: Wow, 78 seqs here! 24 are deleted for not having TM domains. 9 more have no cap helix domain. Sed993Meta19_000000344838_2 and Sed993Meta19_000000439256_2 seen to be truncated before cap helix, but they've got it so we'll include them. Five more have that "kind-of" cap helix domain, which we'll include for now. 44 final.
- 994 Mega: 15 seqs. 6 missing TM domains, 1 missing cap helix. 8 good seqs.
- 994 Meta: 34 seqs. Many look to be junky. Cutting 15 for no TM regions. 5 that have no cap helix domain. 2 of the ones remaining are truncated right before cap helix, but we'll leave them in for now.
- 995 Mega: 4 seqs, all look good.
- 995 Meta: 6 seqs. 3 missing cap helix, one missing TM domains. Only 2 left!
- 996 Mega: 5 seqs. 1 missing cap, one missing TM. Think those two might go together actually, they share 1 reside. Anyways, only 3 left.
- 996 Meta: 10 seqs. 6 are missing TM regions. 1 (Sed996Meta19_000000004556_3) looks like it's one of the unique sequences, one that doesn't align well. This one we'll keep for now. Has 4 TM domains.

I saved a file (`hgcA_good.afa`) with all the keeper sequences, and only kept those with a full cap helix domain and several transmembrane domains.



**Extract gene neighborhood data for hgcA+ scaffolds**

I also wanted to inspect the position and gene neighborhood of the hgcA sequences on their scaffolds.
To do this, I first pulled out the scaffolds for each hgcA+ scaffold, and then the corresponding GFF entries.
I then isolated the gene neighborhoods using a python script I modified from Tyler Barnes.
I pulled out sequences 5000bp from before and after the hgcA sequence, as well as the corresponding GFF entries.
I did this initially for both the raw list of hgcA sequences and saved them here: `hgcA_geneNeighborhood_raw.gff` and `hgcA_geneNeighborhood_raw.fna`.
I also subsetted these two files to get to only include the sequences that passed the trimming (base name `hgcA_geneNeighborhood.fna`).



**Search for downstream hgcB sequences**

First I pulled out the genes downstream from hgcA, by first using a custom python script (`retrieve_downstream_gene_name.py`) to grab the gene names, then pulling out the amino acid sequences.
I searched the corresponding downstream genes using a hgcB HMM I built for the 5M project using the hgcB recovered from my assemblies.
I used 30 as a cutoff score, which I determined from some previous work I had done with it.
I aligned the hits against the HMM, and inspected them in Geneious.
We had 98 total hits.
Of these, 6 of them looked pretty truncated at the C-terminus end, but they did have the compete CM/IECGAC region, so we'll include those.
Six more of them had a partial sequence of that region and four hits had nothing in that region due to being truncated.
Let's look at these in the gene neighborhoods, notes here:  `hgcB_notes.xlsx`.
One of these was truncated at the N-terminus end, but the scaffold wasn't truncated (obviously, since this would mean hgcA is gone).
Think this might be an error with Prodigal, so want to align the DNA segment against the other nucleic acid sequences to see if it does look like hgcB.
I saved out a list of the hgcB names.
For the hgcB.faa file, we'll only include the really good ones, in case someone wants to use these as a seed to search for hgcB.

**Extract depths of hgcA**

Next I wanted to extract the coverage information for all these sequences.
I used the depth function in samtools to calculate the coverage over each nucleotide residue of the hgcA+ scaffolds.
I then used the `calculate_depth_contigs.py` script to calculate the average coverage of each nucleotide over the sequence, excluding 150 bp on both ends of the scaffold, to eliminate any bias against shorter contigs (in case reads aren't mapping to the ends of the scaffold as well).
I did this for each sequence initially identified as hgcA, so that I could see if any of the incomplete sequences are relatively abundant.

I downloaded the depth files to my local computer, then combined them into a single file and normalized them to the single copy gene coverage with the R script `code/assembly_analysis/hgcA/clean_hgcA_coverage.R`.
This saves out a file with the aggregated and normalized depth information here: `dataEdited/assembly_analysis/hgcA/depth/hgcA_coverage_scgNormalization.csv`.



**Clustering of sequences**

Next we'll want to dereplicate the hgcA sequences across our assemblies.
We'll use a cutoff of 97% identity for this process.
In addition to getting the fna file, we'll generate the text file with the cluster information, then download that to my computer.
I also did this using a cut-off of 80% identity to look at patterns of sequence similarity.

**Dereplicate sequences**

Next I wanted to pull all the data I had on hgcA sequences into R to get a good look at them.
I used an Rmd file to keep track of my notes on this: `hgcA_dereplication.Rmd`.
Realized that the 97% cutoff wasn't catching everything here, so dropped it down to 94%.
In the Rmd document, I found that I wasn't always picking the best representative, so I sorted the clusters by length and presence of hgcB and picked a representative that way.
That list got saved out to here: `dataEdited/assembly_analysis/hgcA/dereplication/hgcA_derep_list.txt`.
I uploaded this and used it to pull out the needed hgcA sequences.
Later on, I came back and made a "cluster key" to link all the hgcA seqs in a cluster to the cluster representative.


**Classify hgcA seqs with pplacer workflow**

I also classified the hgcA sequences using Caitlin's pplacer workflow.
The full tutorial can be found here: https://caitlingio.com/tutorial-for-hgcab-amplicon-sequencing-data/
In a nutshell, I align the HgcA sequences to the reference sequences in the refpackage, then use pplacer to place the sequences on the reference tree.
Then we make a sqlite database of the references and use guppy to classify them.
Finally, Caitlin wrote a script that extracts the taxonomic information.
My implementation of this is here: `code/assembly_analysis/hgcA/clean_hgcA_classification.R`.


**Make phylogenetic tree**

*Prepped alignment, preliminary tree*

I generated a consensus alignment between the *hgcA* seqs I found in this study with the Hg-MATE reference database, using MUSCLE.
I then generated a tree using FastTree.
I downloaded the tree file to my local computer, where I visualized it in the `clean_hgcA_tree.R` code file.
There were obviously a ton of unneeded branches.
I manually identified the nodes that needed removal, then used the `tree_subset` function in the `treeio` package to pull out the sequence names for removal.
After consideration, I decided to keep all the isolate sequences in the tree, just for a nice reference.
Read out list here: `dataEdited/assembly_analysis/hgcA/phylogeny/seqs_to_remove.txt`
I used this list to remove the unneeded sequences from my alignment.
To clean the alignment, I masked it in Geneious, trimming out any residue with more than 50% gaps.
The final exported sequence is `hgcA_for_phylogeny_masked.afa`.

*Generate tree in RAxML*

I then generated an ML tree using RAxML.
When doing this, several sequences were found to be identical:
```
IMPORTANT WARNING: Sequences Sed993Mega19_000001657572_1 and Sed993Meta19_000000166988_2 are exactly identical
IMPORTANT WARNING: Sequences Sed993Mega19_000000592489_3 and Sed993Meta19_000000164802_1 are exactly identical
```
There were also many reference sequences that had this warning.
So, I used the reduced file (`hgcA_for_phylogeny_masked.afa.reduced`) for my phylogeny.
Just remember when classifying them to include the two study seqs that were removed.

We used RAxML v8.2.11.
The seed value was 283976.
Rapid bootstrapping was used (with the "-N autoMRE" flag, and a seed number of 2381).
Automatic detection of convergence stopped bootstrapping after 450 replicates.
The phylogeny was generated under a gamma distribution and the best substitution matrix was automatically determined to be LG (-m PROTGAMMAAUTO).

*Clean tree*

The tree was then downloaded and loaded into R (`clean_hgcA_tree.R`).
It was mid-point rooted using the phangorn package and visualized using ggtree.
I replaced the accession numbers of the references with their names, and saved out a cleaned hgcA tree (`dataEdited/assembly_analysis/hgcA/phylogeny/hgcA_clean_tree.rds`).
I'll also use this script to save out a color vector for the tree.
I used these two RDS files, along with the depth information for the hgcA+ scaffolds, to generate an hgcA tree paired with the overall abundance information.


*Generate tree for publications with subset of references using RAxML*

The previous tree that was generated was great for assigning taxonomy, but I wanted an even more cleaned version for publication.
So, I manually added the references to include here: `dataEdited/assembly_analysis/hgcA/phylogeny/final/refs/reference_names_to_use.txt`.
From this, I pulled out the fasta sequences that I needed.
On GLBRC, I aligned the references with the study sequences using MUSCLE.
I masked the alignment at 50% gaps using trimal, then generated a ML tree using RAxML with the usual settings.


**hgcA phylogenetic cluster assignment**

I then assigned a taxonomic group to each of the hgcA sequences, based on the phylogenetic tree that I had built: `results/2019_analysis_assembly/hgcA_tree_RAxML_midpointRooting.pdf`.
This was not the final tree, but was from the original RAxML tree.
That information was saved here: `hgcA_phylogenetic_clusters.xlsx`.


**hgcA abundance analysis**

I looked at the abundance of the hgcA sequences in different ways in the `abundance_hgcA.R` file.


**Bring all hgcA information together**

This is done in `code/assembly_analysis/hgcA/aggregate_hgcA_information.R`.
