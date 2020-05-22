### Protocol for analyzing methanogenesis marker


**Literature notes on methanogenesis analysis**

First I wanted to look into some background reading on how people are identifying methanogenic archaea in metagenomes.

- Boonfei et al: Most of their methanogenesis analysis is in the supplementary info. They compared hydrogenotrophic, acetotrophic, and methylotrophic methanogenesis using PFAM models. MtaB, their marker for methylotrophic methanogenesis was absent. CODH/ACS or AK and PTA are markers for acetotrophic methanogenesis, while the enzymes of the WL pathway in archaea are markers for hydrogenotrophic methanogenesis. Mtr and Mcr are used in both acetotrophic and hydrogenotrophic methanogenesis.

- Carr et al, 2018 assembled a *Methanosaeta* from marine sediments, which was the dominant methanogen. It seemed to be mostly acetoclastic.

I think I'm going to need to search for Mcr, Mtr, and CODH subunits, then generate phylogenies, to see which of them are methanogenic, since all of these subunits can be used in other things.
Let's start with McrA and MtrA.
TIGRFAM has a set of HMMs against different Mtr subunits.
We could use Daan Speth's database from his mcrA paper to build an HMM for that protein.


**mcrA analysis**

I built an HMM for McrA using Daan's database, my notes on it are here: `/Users/benjaminpeterson/Box/references/proteins/methanogenesis/mcrA/mcrA_HMM_construction.md`.
I added this HMM to my metabolic HMM list.
Then, I'll build a tree of mcrA sequences with mine and the sequences from Daan's paper.
In the first iteration of this, the bulk of our sequences clustered closely in the Methanomicrobiales group.
There are only a few references scattered in here, so I'm going to flesh out the reference tree with sequences from refseq as well.
I pulled out 5 that match each sequence, then dereplicated them.
Once I retrieved all those sequences, I concatenated them together and added Daan's sequences and the seqs from this study.
I then generated an alignment of the sequences using muscle.
Finally I generated a tree in FastTree.
I read this into R and visualized the tree.
Looked pretty good, but I had a lot of unneeded sequences.
In R, I read out a list of labels for tips that I didn't need, then pulled those out of the concatenated faa file.

In a new folder, I aligned this subsetted fasta file, then generated a tree using FastTree.
It looked pretty good, so I decided to go ahead with a RAxML run to generate a nicer tree.
First, I masked any residues points with >50% gaps in Geneious, then exported this file.
I uploaded that to GLBRC and generated a tree using RAxML.
