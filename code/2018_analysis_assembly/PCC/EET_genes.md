### Protocol to attempt to find EET genes in assemblies



**Identify potential PCCs**

*Search for MHCs*:
One requisite element of a PCC is the multiheme cytochrome-c (MHCs) protein near the BB-OMP.
To be able to label nearby proteins as MHCs, I searched the entire ORF collection for each assembly using a custom python script developed by Shaomei.
I looked for any protein with 3 or more heme-binding sites.

*Extract proteins to either side*: I then pulled out the ORFs on either side of each MHC.
I then dereplicated this list.

**Search adjacent proteins for BB-OMPs**

Next I searched all these adjacent genes for ones that might be a BB-OMP.
I built my own HMM for this, see notes here: `/Users/benjaminpeterson/Box/references/proteins/EET/EET_reference_notes.md`.
I used this HMM to pull out sequences using a low score threshold (50) then aligned them to the HMM.
They all looked decent.


**Phylogeny**

I wanted to generate a phylogeny to get a better idea of which sequences really look like BB-OMPs.
I generated an alignment (using MUSCLE v3.8.31) of the dereplicated sequences I pulled out with the custom HMM and a cutoff of 50 with the reference PCC BB-OMPs.
I masked the alignment in Geneious at 50% gaps, then generated a tree with FastTree.
I visualized that tree in R: `BBOMP_tree.R`.
The four sequences don't cluster within the BB-OMPs we have as references.
It's likely that they're highly divergent, if they actually are BB-OMPs.


**Extract depth of potential MHCs**

I calculated the depth of coverage over each of these scaffolds, using the same method as for the metabolic proteins.
I cleaned it up using the scripts in `clean_PCC_abundance.R`.
Also in this script is the code to generate a barplot of the abundance.
The coverage never cracks 4 in a sample, so they're pretty low abundance, but we do see a general trend of decreasing coverage of the sequences over time.
I think this is low enough abundance that I won't include it in the metabolic analysis for now, I'll focus on comparing sulfate reduction and methanogenesis.
