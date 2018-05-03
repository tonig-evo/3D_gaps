'''
This programm writes input file for a 1st sites model of Codeml
Input: .paml.fasta
Output: .mod1.paml
'''
import os, sys, re
from Bio import SeqIO
fname = sys.argv[1]
entry_name = re.search(r"(^.+?)\..*", fname).group(1)
outfile = entry_name + '.mod1.paml'
handle = open(fname, "rU") 
seq = list(SeqIO.parse(handle, "fasta")) # read sequences into ordered list
handle.close()

# this peace of code is quite stupid - I tried to catch HS seq from my list
for i in range(len(seq)):
#    print (seq[i].id)
    if seq[i].id == 'Homo_sapiens':
            homo_seq = str(seq[i].seq)

with open(outfile, 'w') as outp:
	first = "  " + str(len(seq)-1) + " " + str(len(homo_seq)) # evaluates # of seq and seq length
#	print first
	outp.writelines(first + '\n')
	for i in range(len(seq)):
                if seq[i].id != 'annotation':
        		outp.write(seq[i].id + '\n')	
	        	outp.write(str(seq[i].seq) + '\n')

