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
In this script, I wanted to root it to match Kartik's tree in his sulfur cycling paper.
First, I read out a PDF of the unrooted tree with the internal nodes labeled (`dataEdited/metagenomes/2019_analysis_assembly/metabolicProteins/sulfur/SRB/dsr/original_dsr_tree.pdf`).
In previous incarnations of looking at this, I've found that if I root by the monophyletic group containing Candidatus_Hydrothermarchaeota_archaeon_JdFR_18_JGI24020J35080_1000005, Aigarchaeota_candidate_division_pSL4_archaeon_ASPF01000004, and Caldivirga_sp.\_MU80 I'll get something that looks like Karthik's tree.
This is node 577.
I can't detect this automatically in R, since even though the tree is read in as unrooted, the MRCA function treats it as though it is rooted, which makes sense.
After rooting the tree with the root function in ape, I read the tree into a PDF (`dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/original_dsr_tree_rooted.pdf`).

I then compared my tree to Karthik's, looking for the rdsr branch.
It looked like I could use the node points "Nitrospirae_bacterium_RBG_19FT_COMBO_42_15" and "RBG_16_scaffold_151951" to encapsulate the reverse dsrA genes (these ones were from bins that were likely oxidizing sulfate).
It was node 712 that isolated this branch.
From these, I subset the tree and isolated the rdsr sequences, then pulled out the 13 remaining dsrA sequences, which we assumed to be the reductive dsrA genes.
I read these out to this list: `dataEdited/2018_analysis_assembly/metabolicProteins/sulfur/dsr/dsrA_red_list.txt`.

I also wanted to built a tree of the reductive dsrA sequences, to get an idea for what taxonomic group these might be associated with and see if there were shifts over the sulfate gradient.
I pulled some references out of NCBI for this, both out of the refseq database and the nr database.
I just midpoint-rooted for this analysis using the phangorn package.

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
