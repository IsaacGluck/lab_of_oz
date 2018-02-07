#!/usr/bin/env python
from itertools import combinations
from ete3 import Tree, TreeStyle
import timeit
from pprint import pprint

# Returns a list of the bootstrap sample trees
# File MUST have a new line at the end
def readTreeFile(filename):
	file = open(filename, "r")
	trees = [] # The list of bootstrap sample trees
	for line in file:
		to_tree = line[:-1]
		trees.append(Tree(to_tree, format=1))
	return trees


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
		temp_list.sort() #ADDED
		to_prune.append(tuple(temp_list))

	return to_prune # a list of tuples (size=4) of leaf names

# TEST combinationsOfLeaves
# t = readTreeFileFirstLine("for_isaac/RAxML_bootstrap.orfg1")
# print len(t)
# print len(combinationsOfLeaves(t))


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

	print "getAllQuartetsAsDictionary started"
	dictonary_of_quartets = {}
	print "File read, dictionary initialized, beginning to loop..."

	combinations_of_leaves = combinationsOfLeaves(trees[0]) # returns list of tuples

	for bootstrap_sample in trees:
		putQuartetsDictionary(bootstrap_sample, combinations_of_leaves, dictonary_of_quartets)

	return dictonary_of_quartets


# TEST getAllQuartetsAsDictionary
# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg3_5")
trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")
# trees = readTreeFile("testTree.txt")
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);


start_time = timeit.default_timer()
dictonary_of_quartets = getAllQuartetsAsDictionary(trees)
end_time = timeit.default_timer()

# # pprint(dictonary_of_quartets)
counter = 0
for key in dictonary_of_quartets:
	print "{0} - {1}: {2}".format(counter, key, dictonary_of_quartets[key])
	counter += 1

print "Dictionary Creation Time: {0}".format(end_time - start_time)