import os, sys
from Bio import SeqIO

indir='Example_Zorro_input/'
alle=os.listdir(indir)
alle.sort()

for a in alle[:]:
	if a+'.zorro' in os.listdir('Gblocks_out'): continue
	os.system('cp '+indir+a+' Gblocks_out/'+a)
	command='./Gblocks Gblocks_out/'+a+' -t=p -k=y -n=y -v=32000 -p=t'
	os.system(command)
