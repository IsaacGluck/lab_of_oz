#!/usr/bin/env python
from itertools import combinations
from ete3 import Tree, TreeStyle
import timeit
from pprint import pprint

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

# Gets every combination of 4 taxa of the tree
def combinationsOfLeaves(tree):
	# Gets all combinations of 4 indices of the tree
	list_of_indices = list(combinations(range(0, len(tree)), 4))
	leaves = tree.get_leaf_names() # List of names of leaves
	to_prune = []

	# converts those indices to actual leaf names
	for inner_list in list_of_indices:
		temp_list = []
		for index in inner_list:
			temp_list.append(leaves[index])
		to_prune.append(tuple(temp_list))

	return to_prune # a list of tuples (size=4) of leaf names

# TEST combinationsOfLeaves
# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")
# print len(t)
# print len(combinationsOfLeaves(t))

# From the list of every combination of 4 taxa (as tuples)
# creates the pruned quartets with matching topologies to the original sample tree
def getQuartets(tree, to_prune):
	list_of_pruned_trees = []

	for list_of_nodes in to_prune:
		t = tree.copy(); # Copy the tree so as not to destroy the parameter
		t.prune(list_of_nodes) # Prune the tree based on the names of the leaves
		list_of_pruned_trees.append(t)

	return list_of_pruned_trees

# TEST getQuartets
# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")
# combos = combinationsOfLeaves(t)
# tree = getQuartets(t, combos)
# print tree




# Returns a list of the bootstrap sample trees
# File MUST have a new line at the end
def readTreeFile(filename):
	file = open(filename, "r")
	trees = [] # The list of bootstrap sample trees
	for line in file:
		to_tree = line[:-1]
		trees.append(Tree(to_tree, format=1))
	return trees

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

# TEST getAllQuartetsAsList
# trees = readTreeFile("for_isaac/RAxML_bootstrap.orfg1")

# c1 = combinationsOfLeaves(trees[0])
# c2 = combinationsOfLeaves(trees[0])

# l1 = getQuartets(trees[0], [c1[0]])
# l2 = getQuartets(trees[0], [c2[0]])

# print l1[0]
# print l2[0]
# print l1[0].compare(l2[0])["rf"]


# test = [trees[0]]
# quartets_for_each_sample = getAllQuartetsAsList(test)
# print quartets_for_each_sample
# print len(quartets_for_each_sample[0])




# Does the heavy lifting for getAllQuartetsAsDictionary
def putQuartetsDictionary(tree, to_prune, dictonary_of_quartets):
	for tuple_of_leaf_names in to_prune:
		t = tree.copy(); # Copy the tree so as not to destroy the parameter
		t.prune(tuple_of_leaf_names) # Prune the tree based on the names of the leaves

		if tuple_of_leaf_names not in dictonary_of_quartets:
			dictonary_of_quartets[tuple_of_leaf_names] = [t, 1, None, 0, None, 0]
		else: # tuple_of_leaf_names in dictonary_of_quartets
			
			# Check 1st Topology
			result0 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][0], unrooted=True)
			if (result0["rf"] == 0): # same topology
				dictonary_of_quartets[tuple_of_leaf_names][1] = dictonary_of_quartets[tuple_of_leaf_names][1] + 1
				continue
			
			# Check 2nd Topology
			if (dictonary_of_quartets[tuple_of_leaf_names][2] == None): # Don't have a 2nd topology yet
				dictonary_of_quartets[tuple_of_leaf_names][2] = t
				dictonary_of_quartets[tuple_of_leaf_names][3] = 1
				continue

			if (dictonary_of_quartets[tuple_of_leaf_names][2] != None): # Already have a 2nd topology
				result1 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][2], unrooted=True) # Check against 2nd topology
				if (result1["rf"] == 0): # same topology
					dictonary_of_quartets[tuple_of_leaf_names][3] = dictonary_of_quartets[tuple_of_leaf_names][3] + 1
					continue

			# Check 3rd Topology
			if (dictonary_of_quartets[tuple_of_leaf_names][4] == None): # Don't have a 3rd topology yet
				dictonary_of_quartets[tuple_of_leaf_names][4] = t
				dictonary_of_quartets[tuple_of_leaf_names][5] = 1
				continue

			if (dictonary_of_quartets[tuple_of_leaf_names][4] != None): # Already have a 3rd topology
				result2 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][4], unrooted=True) # Check against 3rd topology
				if (result2["rf"] == 0): # same topology
					dictonary_of_quartets[tuple_of_leaf_names][5] = dictonary_of_quartets[tuple_of_leaf_names][5] + 1
					continue

			# ERROR
			print "Error: Dictionary is full but topology is not a match"
			print dictonary_of_quartets[tuple_of_leaf_names]
			print t
			break

	return dictonary_of_quartets

# TEST putQuartetsDictionary
# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")
# to_prune = combinationsOfLeaves(t)
# dictonary_of_quartets = {}
# putQuartetsDictionary(t, to_prune, dictonary_of_quartets)
# putQuartetsDictionary(t, to_prune, dictonary_of_quartets)
# pprint(dictonary_of_quartets)
# print len(dictonary_of_quartets)



# Returns a dictionary for each quartet. The key is 
# a tuple of taxa names, value contains the 3 topologies 
# and their counts from the bootstrap sample trees
def getAllQuartetsAsDictionary(trees):
	start_start_time = timeit.default_timer()

	print "getAllQuartetsAsDictionary started"
	dictonary_of_quartets = {}
	print "File read, dictionary initialized, beginning to loop..."

	for bootstrap_sample in trees:
		combinations_of_leaves = combinationsOfLeaves(bootstrap_sample) # returns list of tuples
		putQuartetsDictionary(bootstrap_sample, combinations_of_leaves, dictonary_of_quartets)

	print "Total time: {0}".format(timeit.default_timer() - start_start_time)
	return dictonary_of_quartets


# TEST getAllQuartetsAsDictionary
trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")

# print len(trees[0].get_leaf_names())

dictonary_of_quartets = getAllQuartetsAsDictionary(trees)

dlen = len(dictonary_of_quartets)
total_quartets_analyzed = 0

for key, value in dictonary_of_quartets.iteritems():
	total_quartets_analyzed += (value[1] + value[3] + value[5])


# print total

pprint(dictonary_of_quartets)
counter = 0
for key in dictonary_of_quartets:
	print "{0} - {1}: {2}".format(counter, key, dictonary_of_quartets[key])
	counter += 1