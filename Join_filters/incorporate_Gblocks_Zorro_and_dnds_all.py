import os, sys, random
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid 

# read species to leave out
specs={}
alle=open('excluded_species.txt').readlines()
for a in alle[:]:
	a=a.split(' ')
	if not specs.has_key(a[0]): specs[a[0]]=[]
	specs[a[0]].append(a[1])

indir='../run/'
alle=os.listdir(indir)
alle.sort()

statout=open('statistics.txt','w')

for a in alle[:]:
	#if a in os.listdir('../modified2/'):
	#	#print "hit"
	#	continue
	#if a!='rna100291_P08069_IGF1R': continue
	#if a!='rna101065_P98161_PKD1': continue
	#if a!='rna111162_Q8WY54_PPM1E': continue
	prot = SeqIO.index(indir+a+'/prep/'+a+'.msaprob.out.annot', "fasta")
	cdna = SeqIO.index(indir+a+'/prep/'+a+'.cds.cdna.annot', "fasta")
	newcdna={}
	for p in prot:
		print p, len(prot[p].seq)
		newcdna[p]=''
		curcdna=str(cdna[p[:-2]].seq)
		counter=0
		for let in str(prot[p].seq):
			if let=='-': newcdna[p]+='---'
			else:
				newcdna[p]+=curcdna[counter*3:counter*3+3]
				counter+=1
	
	
		#print p, len(cdna[p[:-2]].seq), len(newcdna[p])
	
	# open Zorro to remove critical positions check with Gblocks
	# Problem, gblocks does not consider all gaps

	gblocks = open('../Gblocks/Gblocks_out/'+a+'-gbMask').read().split('P1;Gblocks')[1]
	#print str(prot['P1;Gblocks'].seq)
	gblocks=gblocks.replace('\n','')
	gblocks=gblocks.replace(' ','')
	
	
	# check for X
	infoX=open(indir+a+'/prep/'+a+'.msaprob.out.annot').read()
	if infoX.find('X')>-1: foundXX=True
	else: foundXX=False
	
	balle=open('../Zorro/Zorro/'+a+'.zorro').readlines()
	critical=''
	print a, len(balle)
	print a, len(gblocks)
	print a, len(str(prot['Homo_sapiens_1'].seq))
	countergap=0
	for bu in range(len(balle)):
		allgap=True
		foundX=False
		
		for p in prot:
			if p=='annotation_1':continue
			testpos=str(prot[p].seq)[bu]
			#if bu==175:
			#	print testpos,p
			if testpos=='X':
				foundX=True
				print "FoundX"
			if testpos!='-':
				allgap=False
				if not foundXX: break
		#if bu< 200: print bu, allgap,str(prot['Homo_sapiens_1'].seq)[bu]
		#if allgap: print "ALLGAP"
		if allgap or foundX:
			countergap+=1
			critical+='0'
			print countergap
		else:
			b=balle[bu]
			if float(b)<9 or gblocks[bu-countergap]!='#':critical+='0'
			else: critical+='1'
	
	#print len(critical)
	#print len(gblocks)
	#print [gblocks]
	#print [critical]
	
	# remove critical positions and sort ordered and disorderd
	disordered=''
	annot=str(prot['annotation_1'].seq)
	for b in range(len(balle)):
		if float(critical[b])==1 and annot[b]=='F' :disordered+='1'
		else: disordered+='0'
	print len(disordered), disordered.count('1')
	dislen=disordered.count('1')
	ordered=''
	for b in range(len(balle)):
		if float(critical[b])==1 and annot[b]=='K' :ordered+='1'
		else: ordered+='0'
	print len(ordered), ordered.count('1')
	ordlen=disordered.count('0')

	# Outwrite files

	# File 1 annotation from alignment
	
	os.system('mkdir ../modified2All/'+a)
	ausgabe=open('../modified2All/'+a+'/'+a+'.annot','w')
	ausgabe.write('>annot\n'+annot+'\n')
	ausgabe.write('>critical\n'+critical+'\n')
	ausgabe.write('>ordered\n'+ordered+'\n')
	ausgabe.write('>disordered\n'+disordered+'\n')
	ausgabe.write('>Homo_sapiens\n'+str(prot['Homo_sapiens_1'].seq)+'\n')
	ausgabe.close()
	
	# File 2 write alignment excluding critical positions and species
	toexclude=[]
	if specs.has_key(a): toexclude=specs[a]
	#print toexclude
	ausgabe=open('../modified2All/'+a+'/'+a+'.prot.ali','w')
	ausgabe2=open('../modified2All/'+a+'/'+a+'.cdna.ali','w')
	specexcluded=0
	specincluded=0
	for p in prot:
		if p[:-2] in toexclude: 
			specexcluded+=1
			continue
			
		if p in ['annotation', 'annotation_1']:
			specexcluded+=1
			continue
		specincluded+=1
		ausgabe.write('>'+p+'\n')
		ausgabe2.write('>'+p+'\n')
		curseq=str(prot[p].seq)
		curdna=newcdna[p]
		seq=''
		seqDNA=''
		for k in range(len(critical)):
			if critical[k]=='1':
				seq+=curseq[k]
				seqDNA+=curdna[k*3:k*3+3]
		ausgabe.write(seq+'\n')
		ausgabe2.write(seqDNA+'\n')
	ausgabe.close()
	ausgabe2.close()
	
	if specincluded==0: continue
	if ordlen==0:continue
	if dislen==0:continue
	# RUN muscle
	os.system('muscle -in ../modified2All/'+a+'/'+a+'.prot.ali -out ../modified2All/'+a+'/'+a+'.prot.muscle')
	
	# Comparison between muscle and msaprot
	align = SeqIO.index('../modified2All/'+a+'/'+a+'.prot.ali' , "fasta")
	muscle = SeqIO.index('../modified2All/'+a+'/'+a+'.prot.muscle', "fasta")
	
	comparison=open('../modified2All/'+a+'/'+a+'.exclude','w')
	specex=0
	allspec=0
	finalincluded=[]
	for m in muscle:
		allspec+=1
		if str(align[m].seq)!=str(muscle[m].seq): 
			comparison.write(m+'\n')
			specex+=1
		else: finalincluded.append(m)
	comparison.close()
	
	# Exclude when homo sapines is missing
	if 'Homo_sapiens_1' not in finalincluded: continue
	statout.close()		
	statout=open('statistics.txt','a')
	statout.write(a+' '+str(allspec-specex)+' '+str(len(critical))+' '+str(ordlen)+' '+ str(dislen)+'\n')
	statout.close()		
	
	# reorder
	finalincluded.remove('Homo_sapiens_1')
	#if len(finalincluded)>29: finalincluded = random.sample(finalincluded, 29)
	finalincluded=['Homo_sapiens_1']+finalincluded
			
	# FINAL WRITING OF PAML FILES
	
	if allspec-specex<10: continue
	
	# disorderd first
	ausgabe2=open('../modified2All/'+a+'/inpaml.All.mod3a','w')
	ausgabe2.write('   '+str(len(finalincluded))+' '+str(disordered.count('1')*3)+'\n')
	for p in finalincluded:
		if p[:-2] in toexclude: continue
		if p in ['annotation', 'annotation_1']:	continue
		if p not in finalincluded: continue
		#specincluded+=1
		#ausgabe.write('>'+p+'\n')
		ausgabe2.write(p[:-2]+'\n')
		curseq=str(prot[p].seq)
		curdna=newcdna[p]
		seq=''
		seqDNA=''
		for k in range(len(critical)):
			if disordered[k]=='1':
				seq+=curseq[k]
				seqDNA+=curdna[k*3:k*3+3]
		#ausgabe.write(seq+'\n')
		ausgabe2.write(seqDNA+'\n')
		
	#ausgabe.close()
	ausgabe2.close()
	# orderd
	ausgabe2=open('../modified2All/'+a+'/inpaml.All.mod3b','w')
	ausgabe2.write('   '+str(len(finalincluded))+' '+str(ordered.count('1')*3)+'\n')

	for p in finalincluded:
		if p[:-2] in toexclude: continue
		if p in ['annotation', 'annotation_1']:	continue
		if p not in finalincluded: continue
		#specincluded+=1
		#ausgabe.write('>'+p+'\n')
		ausgabe2.write(p[:-2]+'\n')
		curseq=str(prot[p].seq)
		curdna=newcdna[p]
		seq=''
		seqDNA=''
		for k in range(len(critical)):
			if ordered[k]=='1':
				seq+=curseq[k]
				seqDNA+=curdna[k*3:k*3+3]
		#ausgabe.write(seq+'\n')
		ausgabe2.write(seqDNA+'\n')
		
	#ausgabe.close()
	ausgabe2.close()
	#all
	ausgabe2=open('../modified2All/'+a+'/inpaml.All.mod1','w')
	ausgabe2.write('   '+str(len(finalincluded))+' '+str(critical.count('1')*3)+'\n')
	for p in finalincluded:
		if p[:-2] in toexclude: continue
		if p in ['annotation', 'annotation_1']:	continue
		if p not in finalincluded: continue
		#specincluded+=1
		#ausgabe.write('>'+p+'\n')
		ausgabe2.write(p[:-2]+'\n')
		curseq=str(prot[p].seq)
		curdna=newcdna[p]
		seq=''
		seqDNA=''
		for k in range(len(critical)):
			if critical[k]=='1':
				seq+=curseq[k]
				seqDNA+=curdna[k*3:k*3+3]
		#ausgabe.write(seq+'\n')
		ausgabe2.write(seqDNA+'\n')
		
	#ausgabe.close()
	ausgabe2.close()
	#print finalincluded
	
	# prune the tree
	specnames=''
	for f in finalincluded:
		specnames+=f[:-2]+' '
	if 'tree.nw_bk' in os.listdir(indir+a+'/prep/'): treefile=indir+a+'/prep/tree.nw_bk'
	else:
		treefile=treefile=indir+a+'/tree.nw'
	#print './nw_prune -v '+treefile+' '+specs
	os.system('./nw_prune -v '+treefile+' '+specnames+' > tree.nw')
	os.system('Rscript R_unroot.R')
	os.system('cp tree_unrooted.txt ../modified2All/'+a+'/tree.30.unrooted.nw')
	
	# create paml inputfiles for model 1 and 3 ( model 2 to be done by Arina)

	
