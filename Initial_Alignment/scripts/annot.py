'''
This script creates fake sequence annotation using information
of disorder locations from MobiDB database as list: [[start, end]]
Input: .msaprob.out and UniProt ID
Output: .msaprob.out.annot .cds.cdna.annot and .annotation (output redirect)
'''
import os, sys, re
from Bio import SeqIO
fname = sys.argv[1]
entry_name = re.search(r"(^.+?)\..*", fname).group(1)
code = sys.argv[2]
cdna_path = '../cdna/'
seq = SeqIO.index(fname, "fasta")
mobid = str(seq[code].seq)
# Create list of disorder ranges from MobiDB disorder annotation
disreg = eval(re.search(r">"+ code+ ".*(\[\[.*\])", open(fname).read(), re.DOTALL).group(1))
print disreg
my_seq = ''
my_DNA = ''
letters = 0
for letter in mobid:
        if (letter == '-'):
                my_seq += '-'
        elif letter.isalpha():
                letters += 1
		in_disreg = 0
                for item in disreg:
                        if letters in range(item[0], item[1]+1):
                                my_seq += 'F'
                                my_DNA += 'TTT'
                                in_disreg = 1
                if not in_disreg:
                        my_seq += 'K'
                        my_DNA += 'AAA'
print 'annotation_1'
print my_seq
print 'annotation'
print my_DNA
# Write new files .msaprob.out.annot and .cds.cdna.annot with fake protein and cdna annotations
with open(fname, 'r') as inpf1, open(entry_name + '.msaprob.out.annot', "w") as outf1:
        outf1.writelines(inpf1)
	outf1.write(">annotation_1\n")
        outf1.write(my_seq + "\n")
with open(cdna_path + entry_name + '.cds', 'r') as inpf2, open(entry_name + '.cds.cdna.annot', "w") as outf2:
	outf2.writelines(inpf2)
        outf2.write(">annotation\n")
        outf2.write(my_DNA + "\n")
# Remove MobiDB sequence from .msaprob.out.annot file
rm_mobidb = "sed -i '/>"+code+"/,/>/{/>/b;d}' " + entry_name + '.msaprob.out.annot'
os.system(rm_mobidb)  # removes MobiDB sequence from prot seq alignment file
rm_mobidb = "sed -i '/>"+code+"/d' " + entry_name + '.msaprob.out.annot'
os.system(rm_mobidb)


