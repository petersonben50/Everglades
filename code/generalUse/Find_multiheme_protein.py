#! /usr/bin/python

'''
This script prints the list of sequences with at least n heme-binding motif "C..CH".
It requires two inputs. The first input is the file you'd like to search, and the
second is the minimum number of heme-binding sites required to include an ORF.
'''

import sys
import re
from Bio import SeqIO

infile = sys.argv[1]

if sys.argv[1] == []:
	print 'Usage:'
	print argv[0] + 'Infile'

nstring = sys.argv[2]
n = int(nstring)
#n = int(input('Please enter the least number of hemes to define "multiheme": '))

heme = re.compile(r"C..CH")

seq_record = list(SeqIO.parse(infile,"fasta"))

i = 0

outfile1 = infile.split('.')[0] + '_' + str(n) + '_heme_count.txt'
outfile2 = infile.split('.')[0] + '_' + str(n) + '_heme_list.txt'

with open(outfile1, 'w') as f:
	f.write('{0}\t{1}\n'.format('gene_oid', 'Num_of_C..CH'))
	for seq in SeqIO.parse(infile,'fasta'):
		heme_sites = re.findall(heme, str(seq.seq))
		if len(heme_sites) >= n:
			i+=1
			f.write('{0}\t{1}\n'.format(seq.id, len(heme_sites)))

with open(outfile2, 'w') as f:
	for seq in SeqIO.parse(infile,'fasta'):
		heme_sites = re.findall(heme, str(seq.seq))
		if len(heme_sites) >= n:
			f.write('{0}\n'.format(seq.id))

print('\nThere are %i sequences and %i of them have at least %i heme-binding sites.' % (len(seq_record), i, n))
print('The results is printed to %s and %s.\n' % (outfile1, outfile2))
