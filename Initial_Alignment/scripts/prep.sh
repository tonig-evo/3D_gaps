#!/bin/bash
#----------------------------------------------------------#
# This is the main script to run file preparation for PAML #
# It takes prot sequence files and cdna sequence files     #
# '.cds' from /prot and /cdna subdirectories, proceeds in  #
# /prep subdirectory and saves output and intermedeate     #
# files in /run subdirectory. Main output files: input seq #
# files for PAML 'mod(1,2,3,3_2).paml', tree.nw and        #
# .annotation                                              #
# Preparation includes: alignment with MSAProb, annotation,#
# pal2nal and tree pruning                                 #
# Run as follows:                                          #
# bash scripts/prep.sh                                     #
#----------------------------------------------------------#

# Paths:
path=$(pwd)
pal2nal='/home/arina/SOFT/pal2nal.v14/pal2nal.pl'
Mammal_tree=$path/Mammal_tree_Toni_names_noroot.tree
 
# 1. Copy prot seq to /prep folder and add MobiDB seq from MobiDB data file 'MobiDB.filtered.fa'
 
mkdir prep
mkdir run

# cd $path/prot
# for file in *.cds
#	do 
#	code=$(echo $file | sed 's/.*_\(.*\)_.*/\1/')
#	cp $file $path/prep/$file.prot
#	grep -A 1 -e $code $path/MobiDB.filtered.fa >> $path/prep/$file.prot
#	done

# 2. Run the procedure of PAML input files preparation for each .cds.prot file in /prep directory

cd $path/prep

sorted_list=($(ls *.cds.prot | sed "s/\.cds\.prot//g")) 
for job in $(seq 0 $((${#sorted_list[@]}-1))); do 
	date1=$(date +"%s")
   	entry=${sorted_list[$job]} 
	echo $entry
	python $path/scripts/align_seq.py $entry.cds.prot 
	# sort prot and cdna seq files for pal2nal
	python $path/scripts/sort_seq.py $entry.msaprob.out.annot $entry.cds.cdna.annot
	$pal2nal $entry.prot.fasta $entry.cdna.fasta -nogap -output fasta > $entry.paml.fasta
	# write PAML input files
	python $path/scripts/write_simple.py $entry.paml.fasta
	python $path/scripts/write_sites.py $entry.paml.fasta
	python $path/scripts/write_dis.py $entry.paml.fasta
	python $path/scripts/write_ord.py $entry.paml.fasta
	# create list to prune the tree
	sp_list=$(cat $entry.paml.fasta | grep -e '^>' | sed 's/^>//' | paste -s -d " " | sed 's/annotation,//')
	$path/scripts/./nw_prune -v $Mammal_tree $sp_list > tree1.nw
	Rscript $path/scripts/unroot.R
	# zip results
	echo $entry.mod1.paml $entry.mod2.paml $entry.mod3.paml $entry.mod3_2.paml $entry.annotation $entry.tree.nw $entry.msaprob.out.annot $entry.cds.cdna.annot $entry.msaprob.out | sed "s/ /\n/g" | zip prep_res -@
	rm $entry.msaprob.out $entry.align.log tree1.nw
	# create subdirectory with the entry name in /run directory and move all fi;es there
	mkdir ../run/$entry
	mkdir ../run/$entry/prep
	mv $entry.annotation ../run/$entry
	mv $entry.mod1.paml ../run/$entry/mod1.paml
	mv $entry.mod2.paml ../run/$entry/mod2.paml
	mv $entry.mod3.paml ../run/$entry/mod3.paml
	mv $entry.mod3_2.paml ../run/$entry/mod3_2.paml
	mv tree.nw ../run/$entry/tree.nw
	mv $entry.msaprob.out.annot ../run/$entry/prep
	mv $entry.cds.cdna.annot ../run/$entry/prep
	mv $entry.prot.fasta ../run/$entry/prep
	mv $entry.cdna.fasta ../run/$entry/prep
	mv $entry.paml.fasta ../run/$entry/prep
	mv $entry.log ../run/$entry/prep
    	wait
	date2=$(date +"%s")
        echo ~~~~~~~~~~~~~~~~~~~~~~~
        echo  Entry $entry was processed for $(($date2-$date1)) s
done


