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





# Reads in a file with a tree structure and returns a ete tree object 
def readTreeFileFirstLine(filename):
	file = open(filename, "r")
	first_tree = file.readline()[:-1] # [:-1] Gets ride of newline at the end of the line
	return Tree(first_tree, format=1)

# PUT IN COMBINATIONSOFLEAVES
# Returns a list of all combinations of the indices of the tree
# def combinationsOfIndices(tree):
# 	return list(combinations(range(0, len(tree)), 4))
# TEST readTreeFileFirstLine and combinationsOfIndices
# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")
# print combinationsOfIndices(t)

# From the list of every combination of 4 taxa (as tuples)
# creates the pruned quartets with matching topologies to the original sample tree
def getQuartets(tree, to_prune):
	list_of_pruned_trees = []

	for list_of_nodes in to_prune:
		t = tree.copy(); # Copy the tree so as not to destroy the parameter
		t.prune(list_of_nodes) # Prune the tree based on the names of the leaves
		list_of_pruned_trees.append(t)

	return list_of_pruned_trees

# TEST getQuartets# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")# combos = combinationsOfLeaves(t)# tree = getQuartets(t, combos)# print tree


# Returns a list for each of the bootstrap samples
# Each list has every combination of 4 taxa
def getAllQuartetsAsList(trees):
	start_start_time = timeit.default_timer()

	quartets_for_each_sample = []
	counter = 0
	for bootstrap_sample in trees:
		start_time = timeit.default_timer()
		
		combinations_of_leaves = combinationsOfLeaves(bootstrap_sample)
		list_of_quartets = getQuartets(bootstrap_sample, combinations_of_leaves)
		quartets_for_each_sample.append(list_of_quartets)
		
		elapsed = timeit.default_timer() - start_time
		print "Loop {0} time: {1}".format(counter, elapsed)
		counter+=1

	print "Total time: {0}".format(timeit.default_timer() - start_start_time)
	return quartets_for_each_sample

# TEST getAllQuartetsAsList# trees = readTreeFile("for_isaac/RAxML_bootstrap.orfg1")# c1 = combinationsOfLeaves(trees[0])# c2 = combinationsOfLeaves(trees[0])# l1 = getQuartets(trees[0], [c1[0]])# l2 = getQuartets(trees[0], [c2[0]])# print l1[0]# print l2[0]# print l1[0].compare(l2[0])["rf"]# test = [trees[0]]# quartets_for_each_sample = getAllQuartetsAsList(test)# print quartets_for_each_sample# print len(quartets_for_each_sample[0])
