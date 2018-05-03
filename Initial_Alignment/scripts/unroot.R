library('ape')
tr <- read.tree("tree1.nw")
unrooted_tr <- unroot(tr)
write.tree(unrooted_tr,"tree.nw")

