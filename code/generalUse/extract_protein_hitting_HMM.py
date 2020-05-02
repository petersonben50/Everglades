
########################################
# extract_proteins_hitting_HMM.py
# Benjamin D. Peterson

# This script pulls out the fasta protein
# sequences that hit the HMM.
########################################

import os
import sys
import Bio
from Bio import SearchIO
from Bio import SeqIO

# Read input arguments from command line into variable names
# Output of the hmmer search
hmmeroutput = sys.argv[1]
# Fasta file from which you need to extract proteins
fasta = sys.argv[2]
# Name of the gene. Really the name of the HMM used.
geneName = sys.argv[3]
# Name of the bin.
binName = sys.argv[4]

# Set name for file to hold scaffold to bin information.
resultsfile = 'temp.' + geneName + '.' + binName + '.fasta'
keyFileName = 'temp.contig_to_bin.' + geneName + '.' + binName + '.tsv'
countFileName = 'temp.count_file.' + geneName + '.' + binName + '.tsv'

# Load up hmmer output
sample = SearchIO.read(hmmeroutput, 'hmmer3-tab')

# Create dictionary of the fasta sequence
faadict = dict()
for seq_record in SeqIO.parse(fasta, "fasta"):
	faadict[seq_record.id] = seq_record.seq


# Set count number
countEm = 0

if sample:
	with open(resultsfile, 'w') as resultFile:
		with open(keyFileName, 'w') as keyFile:
			for sampleID in SearchIO.read(hmmeroutput, 'hmmer3-tab'):
				resultFile.write('>' + str(binName) + '\n' + str(faadict[sampleID.id]) + '\n')
				keyFile.write(str(binName) + '\t' + str(sampleID.id) + '\n')
				countEm = countEm + 1
else:
	print('no hits for ' + fasta)


with open(countFileName, 'w') as countFile:
	countFile.write(binName + '\t' + str(countEm) + '\n')
