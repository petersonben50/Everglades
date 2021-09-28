### Protocol to identify and quantify SCGs in Everglades metagenomes


**Get set up**

First I got the folders set up for the log files from the submission script.
I also set up a list of the SCGs to use.


**Submit jobs**

I wrote these analyses as submission files.
Each job goes after a different gene.
First it pulls the copies of the gene out of all the different assemblies.
It then dereplicates those genes (97% identity).
It then calculates the coverage of all the identified genes in all the different metagenomes.


**Concatenate files and clean up**

Concatenate all the coverage files into one file.
Remove the working directories that I generated.


**Generate normalization vectors**

In the script `code/SCG/SCGs_depth_cleaning.R`, I took the SCG coverage data and used it to generate a normalization vector.
This was done by taking the inverse of the coverage value and multiplying by 100.
That will allow me to multiple coverages by this value to obtain the coverage per 100X SCG coverage.
