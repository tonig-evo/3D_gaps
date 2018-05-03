# Scripts for alignment pipepline used in "Human long intrinsically disordered protein regions are frequent targets of positive selection" (Afanasyeva et al., 2018, Genome research)

## Initial Scripts

The main script to start file preparation - `prep.sh`. Procedure utilises several external programms, such as pal2nal and MSAProb (v. -0.9.7). Before running the script please specify the correct paths in `prep.sh` and `align_seq.py`. Run the pipeline from the parent directory as follows:

  > bash scripts/prep.sh

`prep.sh` evokes following several scripts:
1) `align_seq.py` - main logic of the pipeline is here (alignment, annotation, sequence filtering), this script evokes:
  - MSAProb
  - annot.py
  - final_annot.py
2) sort_seq.py - to sort seq in both protein and cdna files ('Homo_Sapiens' first) for pal2nal procedure
3) pal2nal
4) 4 scripts to write input files for PAML:
  - write_simple.py 
  - write_sites.py
  - write_dis.py
  - write_ord.py
5) nw_prune - to prune a tree according to the list of filtered sequences ('Mammal_tree_Toni_names_noroot.tree' file)

## Added filtering steps
Alignments can be optimised further, python scripts are provided: 
  - realignment (muscle) and comparison with initial alignment
  - pairwise_paml
  - Gblocks (site filter)
  - Zorro (site filter)
