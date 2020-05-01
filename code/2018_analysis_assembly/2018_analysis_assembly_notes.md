#### Notes on overall assembly-based analysis of 2018 samples

**Select assemblies and metagenomes to use**

For the 2018 analysis, I'm going to focus only on the porewater samples, leaving out the surface water samples.
In some preliminary analyses, we didn't find any hgcA sequences in the surface water, which is to be expected since they are mostly oxygenated.
I did not bother assembling surface water samples in this iteration of analysis.

I did assemble the porewater metagenomes with both metaSPADes and MegaHit (individually) and coassembled them with MegaHit.
I included all of the assemblies in my initial analyses, mostly to maximize the net we're casting for hgcA.
This will be especially beneficial for later binning efforts.
The assembly list can be found here: `metadata/lists/2018_analysis_assembly_list.txt`

For the metagenomes to map to these assemblies, I used all the porewater samples.
Surface water samples will be used in binning, but that comes later.
Metagenome list here: `metadata/lists/2018_analysis_assembly_metagenomes_list.txt`
