#!/usr/bin/env python
from ete3 import Tree, PhyloTree
from random import *


# GET THE 1st BOOTSRAP SAMPLE TREE
filename = "for_isaac/RAxML_bootstrap.orfg1"
file = open(filename, "r")
first_tree = file.readline()[:-1] # [:-1] Gets ride of newline at the end of the line

# MAKE IT INTO AN ETE TREE
t = Tree(first_tree, format=1)
print "ORIGINAL TREE\n"
print t

# GET A LIST OF THE LEAVES (by name or node class)
print "\n LEAVES"
# leaves = t.get_leaves()
leaves = t.get_leaf_names()
for index, leaf in enumerate(leaves):
	print (index, leaf)

# GET 4 RANDOM INDICES TO PRUNE
indices = sample(range(0, len(leaves)), 4)
print "\nRANDOM 4 INDICES: " + ', '.join(str(x) for x in indices)

# USE THOSE INDICES TO GET 4 RANDOM NODES
to_prune = []
for index in indices:
	to_prune.append(leaves[index])

print "\nTO PRUNE "
print to_prune
print "\n"

# COPY THE TREE TO NOT LOSE DATA
c = t.copy();

# PRUNE THE TREE
c.prune(to_prune)
print c


# END RESULT
# Old tree still stored in "t"
# Pruned tree stored in "c"