'''
This programm writes input file for a 3d site model of Codeml 
(sites of order only)
Input: .paml.fasta
Output: .mod3_2.paml
'''
import os, sys, re
from Bio import SeqIO
fname = sys.argv[1]
entry_name = re.search(r"(^.+?)\..*", fname).group(1)
outfile = entry_name + '.mod3_2.paml'
handle = open(fname, "rU")
seq = list(SeqIO.parse(handle, "fasta")) # read seq in ordered list
handle.close()

# catch fake annotation sequence
for i in range(len(seq)):
#    print (seq[i].id)
    if seq[i].id == 'annotation':
        annot = str(seq[i].seq)
#print annot

# create list of sequences of ordered sites for all species
code = {}
for i in range(len(seq)):
	code[seq[i].id] = ''
	for j in range(0, len(annot), 3):
		if annot[j:j+3] == 'AAA':
			code[seq[i].id] +=  str(seq[i].seq)[j:j+3]
#	print seq[i].id
#	print code[seq[i].id]
#print len(code), len(code.keys()[0])

with open(outfile, 'w') as outp:
	# evaluate # os seq and seq length for 1st line of file
	first = "  " + str(len(code)-1) + " " + str(len(code.values()[0]))
#	print first
	outp.writelines(first + '\n')
	for i in range(len(seq)): # write sequences of sites of order in right order
                if seq[i].id != 'annotation':
        		outp.write(seq[i].id + '\n')	
	        	outp.write(code[seq[i].id] + '\n')

