import os, sys
from Bio import SeqIO

indir='Example_Zorro_input/'
alle=os.listdir(indir)
alle.sort()

for a in alle[:]:
	command='./zorro_linux_x86_64 -guide trees/'+a+'.tre '+indir+a+' > Zorro/'+a+'.zorro'
	
	# Serial running 
	os.system(command)
	
	# SGE
	#os.system('mkdir scripts')
	#newaus=open('scripts/'+a+'.sh','w')
	#newaus.write(command)
	#newaus.close()
	#os.system('qsub -q standard -e err.e -o out.o -cwd scripts/'+a+'.sh')
	
	
	
