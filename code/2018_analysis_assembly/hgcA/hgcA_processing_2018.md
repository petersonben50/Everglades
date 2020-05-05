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

I saved a file (`hgcA_good.afa`) with all the keeper sequences, and only kept those with a full cap helix domain and several transmembrane domains.



**Extract gene neighborhood data for hgcA+ scaffolds**

I also wanted to inspect the position and gene neighborhood of the hgcA sequences on their scaffolds.
To do this, I first pulled out the scaffolds and the corresponding GFF entries for each hgcA+ scaffold.
I then isolated the gene neighborhoods using a python script I modified from Tyler Barnes.
I pulled out sequences 5000bp from before and after the hgcA sequence, as well as the corresponding GFF entries.
I did this for both the raw list of hgcA sequences as well as the trimmed list.










**Search for hgcA in assembly graphs**

First off, I looked at the assembly graphs of the sequences that were cut off to see if I could lengthen them and get a full length sequence.
There are three such sequences:
- KMBP004E_000000161924_1
    - The corresponding node for this was not attached to any others. Nothing we can do here. hgcA is at the end of the node.
- KMBP004F_000000644814_3
    - The corresponding node for this was not attached to any others. Nothing we can do here. hgcA is at the end of the node.
- KMBP004F_000000253505_6
    - The corresponding node for this was not attached to any others. Nothing we can do here. hgcA is at the end of the node. The hgcA+ node is also scaffolded to another one.

KMBP004F_000000165420 is cut a little short, but it isn't at the end of the node. Might be a mistake in the protein prediction.

All of the scaffolds are on individual nodes except for KMBP004F_000000437519 and KMBP004F_000000037185.
These two are actually in the same portion of the assembly graph.
They both connect to node 1857184499, which is where I'd expect the hgcB sequence to be, since hgcA is at the end.
This is part of NODE_26726_length_13659_cov_9.070130, which appears to be scaffolded to two other contigs, but not connected on the graph.
I pulled out the NA sequences and the GFF files for this scaffold.
In Geneious, I assembled both hgcA+ scaffolds to the downstream one, and the overlap was perfect.
The first gene on the downstream scaffold was hgcB.

I wanted to individually assemble each of these scaffolds.
So, I pulled out the relevant sections of the mapping files and downloaded them to my local computer.
After looking at this, I think it's gonna be tricky to resolve any variations with only a single metagenome that has appreciable coverage of these organisms.
Let's just move on from trying to do this.



**Extract depths of hgcA**

First things first, I wanted to record the depth of each of our metagenomes.
I did this manually by looking through my FastQC files and loading them in `dataEdited/metagenomes/fastQC/metagenome_depth.csv`.

Then I wanted to extract the coverage information for all these sequences.
I used the depth function in samtools to calculate the coverage over each nucleotide residue of the hgcA+ scaffolds.
I then used the `calculate_depth_contigs.py` script to calculate the average coverage of each nucleotide over the sequence, excluding 150 bp on both ends of the scaffold, to eliminate any differences we might see due to being near the end of the scaffold.
I did this for each sequence initially identified as hgcA, so that I could see if any of the incomplete sequences are relatively abundant.

Let's check out the raw coverage data for the incomplete sequences:

- KMBP004E_000000161924_1
    - Abundance of 5.9 in KMBP001E, and 0 in all other metagenomes. Not substantial.
- KMBP004F_000000644814_3
    - Coverage of 2.5 in KMBP004F, essentially 0 in all others.
- KMBP004F_000000253505_6
    - Coverage of 4.8 in KMBP004F, essentially 0 in all others. Also not really worth worrying about.

KMBP004F_000000165420 is also cut a little short, and does count for 9.4X coverage.

I downloaded the depth files to my local computer, then combined them into a single file and normalized them with the R script `clean_hgcA_abundance.R`.




**hgcA dereplication**

Next we'll want to dereplicate the hgcA sequences across our assemblies.
We'll use a cutoff of 97% identity for this process.
In addition to getting the fna file, we'll generate the text file with the cluster information, then download that to my computer.
That file is read in to R (`code/metagenomeAssemblyAnalysis/hgcA/hgcA_percent_identity.R`).

Finally, we'll manually curate the hgcA sequence list that we have.
Most of them should be fine, but there are a few hgcA sequences that cluster with 97% identity, but have less than 95% coverage.
Notes on this process here:
`dataEdited/metagenomes/assemblyAnalysis/hgcA/identification/hgcA_cluster_inspection_notes.ods`.
If a cluster doesn't have >95% cluster coverage, we'll manually look at it.
For these, we'll look at correlation correspondence, the scaffold/gene neighborhood, and the actually sequence identity.
After looking at ours, I think CD-HIT did a fine job of pulling out correct clusters, so we won't change anything.




**hgcA phylogeny**

*Find additional reference sequences*

First we'll want to flesh out our hgcA tree using some reference sequences.
We have our list of 30 hgcA sequences from confirmed methylators, but we'll want to add in a few more.
Let's do that by blasting each of our hgcA sequences against the non-redundant database.
We'll do this using command-line BLAST (v2.9.0), with the blastp function.
The nr database is stored on the GLBRC servers, and was updated on April 7th, 2019.
We'll take the top 6 hits against each sequence, then dereplicate that set of protein sequences against the 30 confirmed methylators.

*Combine reference sequence sets*

We'll then want to combine our reference sequence sets (the confirmed methylators and the blast-derived hgcA sequences), since it is likely that many of these are redundant.
We'll use CD-HIT-2D for this, with a cut-off of 97% identity.
We'll then combine the dereplicated sequences into a single reference file.


*Generate alignment*

First we'll concatenate our sequences with the references that we found.
We'll also add in the hgcA sequences from the Mendota bins.
Once we have our final set of sequences, we'll generate an alignment using MUSCLE (v3.8.31).


*Generate rough tree*

First, to check the dataset, we'll generate a rough tree using FastTree (v2.1.10) to check on the quality of the alignment.
I downloaded this to Geneious and checked it out.
None of the branches look particularly egregiously long, so I went ahead with this full set of sequences.


*Curate alignment*

Before generating a tree in RAxML, I curated the alignment by masking with BMGE1.1 on the [Galaxy portal](https://galaxy.pasteur.fr/?tool_id=toolshed.pasteur.fr%2Frepos%2Fdcorreia%2Fbmge%2Fbmge%2F1.12&version=1.12&__identifer=5anq311b6rv).
I used the BLOSUM30 substitution matrix, with everything else on default (sliding window of 3, a maximum entropy threshold and a gap rate cut-off of 0.5, and a minimum block size of 5).
I downloaded this and saved it as `hgcA_masked.afa`.
Then I checked the alignment in Geneious.
One of the sequences from the assemblies is pretty truncated, as are two references, but I went ahead with the analysis anyways, and will check on those branches afterwards to see how it's affected.
I did no additional trimming on this alignment.


*Generate tree in RAxML*

I then generated an ML tree using RAxML.
We used RAxML (v8.2.11) for this.
Rapid bootstrapping was used (with the "-N autoMRE" flag, and a seed number of 2381).
Automatic detection of convergence stopped bootstrapping after 650 replicates.
The phylogeny was generated under a gamma distribution and the best substitution matrix was automatically determined to be LG (-m PROTGAMMAAUTO).
The tree was then downloaded and loaded into R. It was mid-point rooted using the phangorn package and visualized using ggtree.



**hgcA analysis**

I checked out the phylogeny and the overall depths of the HgcA sequences in `hgcA_phylogeny_abundance.R`.
There were two HgcA sequences that cluster with Geobacterales sequences, one of which was the most abundant sequences we had.
These are the two sequences that are essentially the same, and share a branch on the assembly graph.
We have another sequence that clusters closely with the HgcA sequence from the SYN_0012 bin we recovered in Lake Mendota.
Two more sequences cluster with a group of HgcA seqs from Spirochaetes.
Another sequence is relatively close to the HgcA from *Paludibacter jiangxiensis*.
Finally, we had one HgcA sequence cluster with the HgcA from DES_0019 in Mendota.

In the Excel file `hgcA_groups.xlsx`, I assigned a group name to each of the hgcA sequences, based on the phylogeny of the HgcA sequences.
I used this to generate a stacked bar plot of the *hgcA* coverage in `hgcA_abundance.R`.
The Geobacterales sequences dominate the overall coverage.
