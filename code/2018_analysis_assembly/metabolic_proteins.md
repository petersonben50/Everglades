#### Protocol to identify metabolic proteins in 2019 metagenomes

**Make list of metabolic proteins of interest**

The first thing I did was gather a list of the metabolic proteins I was interested.
This was first developed for the Hells Canyon project, so I copied the files from there.
Files include a csv file with the names and corresponding HMM profiles for each of the metabolic proteins, stored here: `~/Everglades/references/metabolic_HMMs.csv`, with the HMMs stored in the corresponding folder (`~/Everglades/references/metabolic_HMMs`).


**Identify potential metabolic genes**

For each sequence of interest, I searched each assembly with the corresponding HMM using hmmsearch in hmmer (v3.2.1).
I used the trusted cut-off score for each HMM.
I extracted the protein sequences, then concatenated all the sequences from all the assemblies and aligned them to the HMM of interest using hmmalign in hmmer.


**Extract depths of all scaffolds**

First, I extracted the depths of all scaffolds that might have had a metabolic gene of interest on them.
I concatenated the unique scaffold names to a file, then use that file to pull out the depth of all these scaffolds of interest in all the metagenomes.
I'll pull out the depths of the scaffolds of confirmed sequences later.
To do this, I used the depth function on samtools (v1.10) to get the depth at each residue.
I then used a custom python script to calculate the average coverage across each scaffold of interest within a metagenome.
I didn't include the coverages of the residues within 150 residues of the end of the scaffold, just to reduce any noise that might come from being near the end of the scaffold.
One thing that this doesn't account for is scaffolding gaps, where the coverage tends to drop dramatically... Might need to think about that in the future.
Finally, I downloaded the depth files to my local computer and aggregated the depths across all metagenomes using R: `code/2018_analysis_assembly/metabolicProteins/metabolic_proteins_depth_aggregate.R`.


**Dereplicate genes**

I then used CD-HIT (v4.6) to dereplicate the metabolic genes, leaving us with a final set of sequences to work with.
These sequences will then also be aligned to their respective HMM.
These will be the files I use for further analysis.
