
#### Assemble needed sequences using EMIRGE


**Set up**

First, I need to set up a file with the mean insert sizes and the corresponding SD values.
I pulled the mean insert size data from the library prep sheet here: `protocols/metagenomeSequencing/2018_trip/20180921_SeqProduction_McMahon.xlsx`.
I saved this into a tab-separated file.
We'll save that here: `dataEdited/16S_from_metagenomes_2018/setup/MG_insert_sizes_PW.txt`.
Unfortunately we don't have standard deviations for these values.
In our HCC data, there was an average of 27% coefficient of variation among the insert sizes.
We'll use that average for this, and hope for the best!
Might be worth doing a few iterations of this.


**Run EMIRGE on metagenomes**

I individually ran EMIRGE on each of the porewater MGs from 2018, which is stored here: `~/Everglades/metadata/lists/2018_analysis_assembly_metagenomes_list.txt`.
I unzipped the fastq files temporarily for this analysis, then deleted them once it was done.
I only used the forward and reverse reads in this analysis.
The analysis was based on the SILVA database (SILVA_138_SSURef_NR99_tax_silva_trunc).
The insert size and standard deviation was called from the tab-separated sheet that I had previously set up (see above).





**Pull out all FASTA sequences**

Once EMIRGE was done, I pulled out all the fasta sequences that it had assembled and renamed them using the `emirge_rename_fasta.py` python script in EMIRGE.
This add the prior abundances to the names.
I then added the metagenome ID to each name to ensure they were unique and moved them to a single folder.


**Cluster 16S sequences across metagenomes**

I only wanted one set of 16S sequences, so I used CD-HIT to cluster the fasta files.
I used a cut-off of 0.97 to cluster them, then generated a key file, where the 16S sequence name and the representative sequences are stored in two separate columns.


**Extract abundance of 16S sequences**

I based this section on the workflow from [Vincent Denef's group](https://github.com/DenefLab/EMIRGE).
I extracted the total number of reads mapped to each metagenome's set of 16S sequences using samtools.
I also pulled out the fractional abundance of each sequence from the name of the 16S sequence (stored as "NormPrior"), and then cleaned up the count file.


**Rename seqs in abundance file with representative ID**

I then took the relative abundance file I generated in the previous step and renamed each sequence with the representative sequence generated using CD-HIT.


#### Classify 16S sequences using mothur

I then classified all of my 16S sequences using [TaxAss](https://github.com/McMahonLab/TaxAss).
This is also outlined on the Denef github workflow, cited above.
I used the FreshTrain dataset for my primary database (accessed March 10th, 2019) and Silva for the secondary one (also accessed March 10th, 2019).
I won't outline all the details here since I've mostly just followed the usual workflow for this.
