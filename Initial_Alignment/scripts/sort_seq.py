'''
This programm sorts both protein and cDNA sequences files
in same order, 'Homo_Sapiens' sequence first.
Input files: .msaprob.out.annot .cds.cdna.annot
Output files: .prot.fasta .cdna.fasta
'''
import sys, re
from Bio import SeqIO
protf = sys.argv[1]
seqf = sys.argv[2]
entry_name = re.search(r"(^.+?)\..*", protf).group(1)

prot = SeqIO.index(protf, "fasta")
cdna = SeqIO.index(seqf, "fasta") 
out_prot=open(entry_name+'.prot.fasta','w')
out_cdna=open(entry_name+'.cdna.fasta','w')

out_prot.write('>Homo_sapiens\n'+str(prot['Homo_sapiens_1'].seq)+'\n')
out_cdna.write('>Homo_sapiens\n'+str(cdna['Homo_sapiens'].seq)+'\n')	


for p in prot:
	if p=='Homo_sapiens_1':continue
	out_prot.write('>'+p[:-2]+'\n'+str(prot[p].seq)+'\n')
	out_cdna.write('>'+p[:-2]+'\n'+str(cdna[p[:-2]].seq)+'\n')	

out_prot.close()
out_cdna.close()

