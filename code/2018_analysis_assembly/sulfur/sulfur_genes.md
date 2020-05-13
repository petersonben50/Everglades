#### Protocol to confirm sulfur cycling genes in metagenome assembly analysis

**Sulfate reducers**

I checked the abundance of the genes of interest in the R script `code/2019_analysis_assembly/sulfur/SRB_abundance.R`
For SRBs, I looked at the dsr complex (A and D subunits), and the anaerobic sulfite reductase (asr) complex.

*dsrA*

First I wanted to look at the dsr genes I pulled out.
There are many reverse dsr genes (rdsr) for sulfide oxidation, so I needed to separate those out from the reductive ones.
I looked at the dsrA gene alignment in Geneious, they all looked good, so going forward with them.
I then aligned that alignment to Karthik's dsrA gene set from his sulfur paper.
He sent me his set of reference sequences on March 23rd, which I cleaned up (`/Users/benjaminpeterson/Box/references/proteins/dsr/dsrA/dsrA_karthik_notes.md`).
I masked the alignment on 50% gaps in Geneious.
I then generated a ML tree using RAxML (v8.2.11), with a gamma distribution and autodetection of the protein substitution matrix.
I downloaded `RAxML_bipartitions.dsrA` to the local computer (`dataEdited/metagenomes/2019_analysis_assembly/metabolicProteins/sulfur/SRB/dsr`).

I started looking at the tree here: `code/2019_analysis_assembly/sulfur/dsr_tree.R`.
In this script, I read in the tree, and midpoint rooted it.
I read out the tree to `dataEdited/metagenomes/2019_analysis_assembly/metabolicProteins/sulfur/SRB/dsr/original_dsr_tree.pdf` and manually inspected the tree.
It looked like I could use the node points "Nitrospirae_bacterium_RBG_19FT_COMBO_42_15" and "rifcsplowo2_02_scaffold_46312" to encapsulate the reverse dsrA genes (these ones were from bins that were likely oxidizing sulfate).
From these, I subset the tree and isolated the rdsr sequences.
Looks like I only had 3 reductive dsrA genes.
I read these out to this list: `dataEdited/metagenomes/2019_analysis_assembly/metabolicProteins/sulfur/SRB/dsr/dsrA_red_list.txt`.

I then read this list into `code/2019_analysis_assembly/sulfur/SRB_abundance.R` and found the corresponding depths for each of the scaffolds.
I plotted out the depth of dsrA in each sample.
While it's mostly present at RM300, one of the dsrA sequences (some Deltaproteobacteria) has some coverage at 305/310.



*dsrD*

To double-check our abundance of sulfate-reducers, I also looked at the abundance of dsrD.
For this, I downloaded the `dsrD_derep_list.txt` to my local computer.
Clear gradient in dsrD levels from high sulfate to low sulfate samples.


*asr*

There are no anaerobic sulfite reductases.



**Sulfur oxidizers**

Next we'll look for sulfur oxidizers.
I used two metrics: looking for the sox genes, and looking at the rdsr genes.

*sox genes*

Looked at the alignments in Geneious first.
They looked to be good, but the numbers are very inconsistent.
A and X have over 300 hits, while B and C have around 120 hits.
A and X share 4 hits, but no shared hits between other genes.
Not sure what to make of this.
I just used all of the hits that we got from this group.

*rdsr genes*

I also wanted to look at the abundance for reverse dsrA as a proxy for sulfide oxidation.
I got this list out of the dsr_tree.R script.

*Overview of SOBs*

Nothing too definitive here.
RM286 has the highest coverage of sulfur oxidizing organisms, but it then increases from RM300 up to RM318.
Not a clear answer here.
