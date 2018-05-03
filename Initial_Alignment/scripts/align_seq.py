#---------------------------------------------------------#
# This script contains main logic of PAML preparation     #
# pipeline: alignment, annotaion and sequences filtering  #
# Input file: .cds.prot                                   #
# Output files: .msaprob.out.annot, .cds.cdna.annot,      #
#  .annotation						  #
#---------------------------------------------------------#

import os, sys, re, subprocess, shlex, operator
from Bio import SeqIO
fname = sys.argv[1]
entry_name = re.search(r"(^.+?)\..*", fname).group(1)
sys.stdout = open(entry_name+".log", "w")
code = re.search(r"_(.*)_", fname).group(1)

my_path = '/home/arina/test_run'
msa_path = '/home/arina/SOFT/MSAProbs-0.9.7/MSAProbs/./msaprobs'
annot_path = my_path + '/scripts/annot.py'
final_annot_path = my_path + '/scripts/final_annot.py'

print "Entry: ", entry_name

def run_alignment():
#
# Evokes MSAProb for sequence alignment 
# input file: .cds.prot (contains MobiDB sequence)
# output file: .msaprob.out in fasta format
#
	align = msa_path + ' ' + fname + ' -v -o ' + entry_name + '.msaprob.out'
	command = shlex.split(align)
	with open(entry_name+'.align.log', 'w') as logfile:
		p = subprocess.Popen(command, stderr=logfile)	
		print 'Alignment in progress...'
		p.communicate()

def run_annot():
#
# Evokes annotation procedure annot.py 
# input: .msaprob.out and UniProt ID
# output: .msaprob.out.annot (contains fake annotation, MobiDB sequence is removed)
#         .cds.cdna.annot (cdna sequences + fake annotation)
#
	annotate = 'python ' + annot_path + ' ' + entry_name + '.msaprob.out' + ' ' + code
	command = shlex.split(annotate)
	with open(entry_name+'.annotation', 'w') as logfile2:
		p = subprocess.Popen(command, stdout=logfile2)  
		print 'Annotation in progress...'
		p.communicate()

def stats():
#
# Calculates statistics for each of the sequences for further seq filtering:
# stat1: contentof Homo Sapiens seq sites in sequences for each species
# stat2: content of disorder sites in sequences for each species
# stat3: content of HS sites in final (after pal2nal) HS sequence
# stat3_dis: content of disorder sites in final HS sequence
#
	global stat1, stat2, stat3, stat3_dis, all_homo, dis_homo
	prot = SeqIO.index(entry_name + '.msaprob.out.annot', "fasta")
        homo_seq = str(prot['Homo_sapiens_1'].seq)
        annot = str(prot['annotation_1'].seq)
        dis_homo = 0
	all_homo = sum(c != '-' for c in homo_seq)
	for i in range(0, len(homo_seq)):
        	if homo_seq[i].isalpha():
#	               	print i, homo_seq[i]
	                if annot[i] == 'F':
        	                dis_homo += 1
	stat1 = dict.fromkeys(prot.keys(), 0)
	stat2 = dict.fromkeys(prot.keys(), 0)
	stat3 = 0
	stat3_dis = 0
	for i in range(0, len(homo_seq)):
		if homo_seq[i].isalpha():
			gap = 0
			for p in prot:
				if str(prot[p].seq)[i] == '-':
					gap = 1
					break
			if not gap:
				stat3 += 1
				if annot[i] == 'F':
					stat3_dis += 1
	for p in prot:
		length = sum(c != '-' for c in str(prot[p].seq))
	        for i in range(0, len(homo_seq)):
        	        if homo_seq[i].isalpha():
                	        if prot[p].seq[i].isalpha():
                        	        stat1[p] += 1
                                	if annot[i] == 'F':
                                        	stat2[p] += 1

def insertion():
#
# Check for long insertions in sequences: if more than
# 20% of sites don't have homologous sites in HS seq - exclude sequence
# 
	global change_check
	prots = SeqIO.index(entry_name + '.msaprob.out.annot', "fasta")
	homo_seq = str(prots['Homo_sapiens_1'].seq)
	for p in prots:
		insert = 0
		prot = str(prots[p].seq)
		prot_len = sum(c != '-' for c in prot)
		for i in range(len(prot)):
			if (prot[i].isalpha() and homo_seq[i] == '-'):
				insert += 1
#		print p, ': ', insert/float(prot_len)*100, '%'
		# Remove bad sequence from .cds.prot and .msaprob.out.annot files
		if insert/float(prot_len)*100 > 20.0:
			if p != 'annotation_1':
				print "Sequence with big insertion: ", p
				rm_seq = "sed -i '/>"+p+"/,/>/{/>/b;d}' " + fname
                	        os.system(rm_seq)  # removes sequence from prot seq alignment file
				rm_seq = "sed -i '/>"+p+"/d' " + fname
                	        os.system(rm_seq)
				rm_seq = "sed -i '/>"+p+"/,/>/{/>/b;d}' " + entry_name + '.msaprob.out.annot'
        	                os.system(rm_seq)  # removes sequence from prot seq alignment file
	                        rm_seq = "sed -i '/>"+p+"/d' " + entry_name + '.msaprob.out.annot'
                	        os.system(rm_seq)
				change_check = 1

def mismatch():
#
# Check for mismatches in HS and MobiDB sequences: places gaps in mismatched 
# sites in HS seq
#
	global mism
        prot2 = SeqIO.index(entry_name + '.msaprob.out', "fasta")
        homo_seq = str(prot2['Homo_sapiens_1'].seq)
        mobid = str(prot2[code].seq)
	cdna = SeqIO.index(entry_name + '.cds.cdna.annot', "fasta")
	homo_cdna = str(cdna['Homo_sapiens'].seq)
	new_hs = ''
	new_cdna = ''
#        print 'Mismatch length: ', len(homo_seq), len(mobid)
	letters = 0
	mism = 0
	for i in range(0, len(homo_seq)):
                if homo_seq[i] == mobid[i]:
			new_hs += homo_seq[i]
			if homo_seq[i].isalpha():
				new_cdna += homo_cdna[3*letters:3*letters+3]
				letters += 1
		else:
			new_hs += '-'
			mism += 1
			if homo_seq[i].isalpha():
				letters += 1
#	print new_hs == homo_seq
#	print 'Old HS: ', homo_seq
#	print 'New HS: ', new_hs
#	print 'Old cDNA: ', homo_cdna
#	print 'New cDNA: ', new_cdna
	print "Number of mismatches in HS amd MobiDB: ", mism
	# Write new HS prot and cdna sequences in .msaprob.out.annot and .cds.cdna.annot files
	if new_hs != homo_seq:	
		rm_hs = "sed -i '/>Homo_sapiens_1/,/>/{/>/b;d}' " + entry_name + '.msaprob.out.annot'
		os.system(rm_hs)	# removes old HS sequence from prot seq alignment file
		rm_hs = "sed -i '/>Homo_sapiens_1/d' " + entry_name + '.msaprob.out.annot'
		os.system(rm_hs)
		with open(entry_name+'.msaprob.out.annot', 'a') as outf1:
		        outf1.write(">Homo_sapiens_1\n")
		        outf1.write(new_hs + "\n")
                rm_hs = "sed -i '/>Homo_sapiens/,/>/{/>/b;d}' " + entry_name + '.cds.cdna.annot'
                os.system(rm_hs) 	# removes old HS cDNA sequence from cDNA seq alignment file
                rm_hs = "sed -i '/>Homo_sapiens/d' " + entry_name + '.cds.cdna.annot'
                os.system(rm_hs)
                with open(entry_name+'.cds.cdna.annot', 'a') as outf2:
                        outf2.write(">Homo_sapiens\n")
                        outf2.write(new_cdna + "\n")


STAT3 = 0
STAT3_dis = 0
# 1. Run 1st alignment
run_alignment()
# 2. Run 1st annotation procedure
run_annot()
change_check = 0
# 3. Calculate statistics and filter sequences
while STAT3 < 80 or STAT3_dis < 80:
	stats()	
	stat_sort = []
	for p in stat1.keys():
		item = (p, stat1[p], stat2[p])
		stat_sort.append(item) # list of distionaries with two statistics 
#	print stat_sort
	for item in sorted(stat_sort, key = lambda x: (x[2], x[1]), reverse=True): # sort seq first by dis then by HS content
		p = item[0]
        	STAT1 = stat1[p]/float(all_homo)*100       # statistics of HS sites in each seq
		try:
			STAT2 = stat2[p]/float(dis_homo)*100       # statistics of disorder sites in each seq
		except ZeroDivisionError:
			break
        	print p, "\nHS: %6.2f,\t dis: %6.2f" % (STAT1, STAT2)
	STAT3 = stat3/float(all_homo)*100       # statistics of HS sites after gap removal
	try:
		STAT3_dis = stat3_dis/float(dis_homo)*100       # statistics of disorder sites after gap removal
	except ZeroDivisionError:
		break
	print "HS_final: %6.2f,\t dis_final: %6.2f" % (STAT3, STAT3_dis)
	# Remove sequence with bad statistics
        if STAT3 < 80 or STAT3_dis < 80:
        	to_remove = sorted(stat_sort, key = lambda x: (x[2], x[1]), reverse=True)[-1][0]
	        print "~~~~~~~~~~\nSeq to remove: ", to_remove
		if to_remove != 'annotation_1':
	        	rm_seq = "sed -i '/>"+to_remove+"/,/>/{/>/b;d}' " + fname
        		os.system(rm_seq)  # removes sequence from prot seq alignment file
        		rm_seq = "sed -i '/>"+to_remove+"/d' " + fname
        		os.system(rm_seq)
			rm_seq = "sed -i '/>"+to_remove+"/,/>/{/>/b;d}' " + entry_name + '.msaprob.out.annot'
                        os.system(rm_seq)  # removes sequence from prot seq alignment file
                        rm_seq = "sed -i '/>"+to_remove+"/d' " + entry_name + '.msaprob.out.annot'
                        os.system(rm_seq)
			change_check = 1
		else:
			break
# 4. Check for long insertions and filter sequences
insertion()
# 5. If some sequences were removed repeat alignment and annotation
if change_check == 1:
	run_alignment()
	run_annot()
# 6. Check for mismatches of HS and MobiDB sequences and rewrite HS seq in case 
# of mismatches with gaps in mismatched sites
mismatch()
# 7. Create final fake annotation with gaps 'final_annot.py'
final_annot = "python "+final_annot_path+" "+entry_name 
os.system(final_annot)
# 8. Calculate final statistics if HS seq was changed due to mismatches with MobiDB seq
if mism > 0:
	stats()
        STAT3 = stat3/float(all_homo)*100       # statistics of HS sites after gap removing
	try:
		STAT3_dis = stat3_dis/float(dis_homo)*100       # statistics of disorder sites after gap removing
	except ZeroDivisionError:
		print "Error: Homo_Sapiens sequence doesn't contain disorder sites!!!"
print "~~~~~~~~~\nFinal statistics: HS %6.2f,\t dis %6.2f" % (STAT3, STAT3_dis)

