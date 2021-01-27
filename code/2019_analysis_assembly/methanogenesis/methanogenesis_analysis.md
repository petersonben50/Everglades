### Protocol for analyzing methanogenesis marker


**Literature notes on methanogenesis analysis**

First I wanted to look into some background reading on how people are identifying methanogenic archaea in metagenomes.
See notes on this in 2018 assembly analysis.

**mcrA analysis**

Using a personal mcrA HMM: `/Users/benjaminpeterson/Box/references/proteins/methanogenesis/mcrA/mcrA_HMM_construction.md`.
I pulled out five references sequences from NCBI RefSeq database that match each sequence, then dereplicated them.
I also dereplicated them against Daan's sequences.
Once I retrieved all those sequences, I concatenated them together and added Daan's sequences and the seqs from this study.
I then generated an alignment of the sequences using muscle.
I downloaded this alignment and masked it at 50% gaps.
Finally I generated a tree in FastTree.
I read this into R and visualized the tree: `code/2019_analysis_assembly/methanogenesis/methanogen_trees.R`.
Looks pretty good, it could be cleaned up but doesn't really need to be for our purposes.

The *mcrA* sequences fell into three general clusters.
I compared the tree to the tree in Daan's paper and classified the three clusters based on his classification.
This data can be found here: `dataEdited/2019_analysis_assembly/metabolicProteins/methanogenesis/mcrA/phylogeny/mcrA_clusters.xlsx`.

In R, I first looked at the overall coverage of *mcrA*:
`code/2019_analysis_assembly/methanogenesis/methanogen_abundance.R`
Then I looked at the abundance while considering the taxonomic composition.

**methanogenic marker abundance**

I also looked at the abundance of 8 different markers of methanogenesis.
These markers were identified by sequence homology using confirmed methanogens and have no function associated with them.
I looked at the abundance in the same R file.
These show the same abundance pattern as *mcrA*.
