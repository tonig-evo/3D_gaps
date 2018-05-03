'''
This script creates final fake annotation after all procedures,
immitates pal2nal work to place gaps in annotation
'''
import os, sys
from Bio import SeqIO

entry_name = sys.argv[1]

file_cdna=entry_name+'.cds.cdna.annot'
file_prot=entry_name+'.msaprob.out.annot'

cdnas = SeqIO.index(file_cdna, "fasta")
prots = SeqIO.index(file_prot, "fasta")
stopp=["TAA","TAG","TGA"]

# Create list of sites to exclude from annotation: gaps and stop 
# codon positions
exclude=[]
for r2 in prots:
	dna=str(cdnas[r2[:-2]].seq)
	prot=str(prots[r2].seq)
	gaps=0
	for r in range(len(prot)):
		if prot[r]=='-': 
			gaps+=3
			exclude.append(r)
		elif dna[3*r-gaps:3*r+3-gaps] in stopp: exclude.append(r)

exclude=list(set(exclude))
exclude.sort()
# print exclude

new_hs=''
new_annot=''
annot=str(prots['annotation_1'].seq)
homo=str(prots['Homo_sapiens_1'].seq)
homo_dna=str(cdnas['Homo_sapiens'].seq)
annot_dna=str(cdnas['annotation'].seq)
gaps=0
for r in range(len(homo)):
        if homo[r]=='-':
                gaps+=3
                new_hs+='---'
                #exclude.append(r)
        elif r not in exclude:
                new_hs+=homo_dna[3*r-gaps:3*r+3-gaps]
        elif r in exclude:
                new_hs+='---'
gaps=0
for r in range(len(annot)):
        if annot[r]=='-':
                gaps+=3
                new_annot+='---'
                #exclude.append(r)
        elif r not in exclude:
                new_annot+=annot_dna[3*r-gaps:3*r+3-gaps]
        elif r in exclude:
                new_annot+='---'


#print new_hs
#print new_annot

# Append new annotation in .annotation file
with open(entry_name+'.annotation', 'a') as annotfile:
	annotfile.write('PAML.annot\n')
        annotfile.write(new_annot+'\n')

