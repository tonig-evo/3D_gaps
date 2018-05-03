import os, sys

indir=# Your directory
alle=os.listdir(indir)

alle.sort()
for a in alle[:]:
	if a+'.pair' in os.listdir('results'): continue
	# This copies the disordered region under consideration
	if 'codeml_in' not in files:
		os.system('cp '+indir+a+'/mod3.paml codeml_in')
	os.system('./codeml codeml_pairwise.ctl')
	os.system('mv codeml_pair results/'+a+'.pair')
	os.system('rm codeml_in')


