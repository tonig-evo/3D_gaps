'''
This programm writes input file for a 2d site specific model of Codeml
It creates site annotation with numbers where 1 indicates sites of
order and 2 - sites of disorder
Input: .paml.fasta
Output: .mod2.paml
'''
import os, sys, re
from Bio import SeqIO
fname = sys.argv[1]
entry_name = re.search(r"(^.+?)\..*", fname).group(1)
outfile = entry_name + '.mod2.paml'
outfile2 = entry_name + '.annotation'
handle = open(fname, "rU")
seq = list(SeqIO.parse(handle, "fasta")) # read sequences in ordered list
handle.close()

# find annotation
for i in range(len(seq)):
    if seq[i].id == 'annotation':
            annot = str(seq[i].seq)
code = ''

# create site annotation with 1's for sites of order and 2's for disorder
for i in range(0, len(annot), 3):
#	print annot[i:i+3]
	if annot[i:i+3] == 'AAA':
		code += '1'
	elif annot[i:i+3] == 'TTT':
		code += '2'

with open(fname, 'r') as inp, open(outfile, 'w') as outp:
    first = "  " + str(len(seq)-1) + "  " + str(len(annot)) + "  G" # evaluate # of seq and seq length 
#    print first
    outp.writelines(first + '\n') # write 1st line
    outp.writelines('G 2\n') # write 2d line 
    outp.write(code + '\n') # write site annotation
    for i in range(len(seq)): # write sequences in right order
        if seq[i].id != 'annotation':
            outp.write(seq[i].id + '\n')
            outp.write(str(seq[i].seq) + '\n')
