#!/usr/bin/env python
from itertools import combinations
from ete3 import Tree, TreeStyle
from math import log
from pprint import pprint
import timeit
import dendropy

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
# 	return to_prune # a list of tuples (size=4) of leaf names
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

	dictonary_of_quartets = {}

	# Adds all topologies to the dictionary with a default value of 0
	for tuple_of_leaf_names in to_prune:
		t1 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[1], tuple_of_leaf_names[2], tuple_of_leaf_names[3])
		t2 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[2], tuple_of_leaf_names[1], tuple_of_leaf_names[3])
		t3 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[3], tuple_of_leaf_names[1], tuple_of_leaf_names[2])
		dictonary_of_quartets[tuple_of_leaf_names] = [Tree(t1), 0, Tree(t2), 0, Tree(t3), 0]

	return dictonary_of_quartets # a list of tuples (size=4) of leaf names


# Does the heavy lifting for getAllQuartetsAsDictionary
def putQuartetsDictionary(tree, dictonary_of_quartets):
	for tuple_of_leaf_names in dictonary_of_quartets: # tuple_of_leaf_names is the key
		t = tree.copy(); # Copy the tree so as not to destroy the parameter
		t.prune(tuple_of_leaf_names) # Prune the tree based on the names of the leaves


		# Check 1st Topology
		result0 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][0], unrooted=True)
		if (result0["rf"] == 0): # same topology
			dictonary_of_quartets[tuple_of_leaf_names][1] = dictonary_of_quartets[tuple_of_leaf_names][1] + 1
			continue

		# Check 2nd Topology
		result1 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][2], unrooted=True) # Check against 2nd topology
		if (result1["rf"] == 0): # same topology
			dictonary_of_quartets[tuple_of_leaf_names][3] = dictonary_of_quartets[tuple_of_leaf_names][3] + 1
			continue

		# Check 3rd Topology
		result2 = t.compare(dictonary_of_quartets[tuple_of_leaf_names][4], unrooted=True) # Check against 3rd topology
		if (result2["rf"] == 0): # same topology
			dictonary_of_quartets[tuple_of_leaf_names][5] = dictonary_of_quartets[tuple_of_leaf_names][5] + 1
			continue

		# ERROR
		print "Error: Topology is not a match"
		print dictonary_of_quartets[tuple_of_leaf_names]
		print t
		break

	return dictonary_of_quartets


# Returns a dictionary of each quartet. The key is
# a tuple of taxa names, value contains the 3 topologies
# and their counts from the bootstrap sample trees
def getAllQuartetsAsDictionary(trees):

	print "getAllQuartetsAsDictionary started"
	dictonary_of_quartets = combinationsOfLeaves(trees[0]) # returns dictonary_of_quartets
	print "File read, dictionary initialized, beginning to loop..."

	for bootstrap_sample in trees:
		putQuartetsDictionary(bootstrap_sample, dictonary_of_quartets)

	return dictonary_of_quartets


# Takes in a list of dictionaries from M gene trees, merges them and normalizes values
# Also calculates the Internode Certainty Value (based on Shannon's Entropy)
# Returned dictonary with tuples of quartets as the keys and lists as the values -
# The lists are structured as follows [t1, P(t1), t2, P(t2), t3, P(t3), IQ]
def mergeQuartetDictionaries(list_of_dictionaries, bootstrap_cutoff):
	print "Merging {0} dictionaries...".format(len(list_of_dictionaries))

	final_dictionary = {}

	while len(list_of_dictionaries) != 0:
		current_dict = list_of_dictionaries.pop(0)

		if (len(current_dict) == 0): # make sure dictionary is not empty
			continue

		for key in current_dict:
			value = current_dict[key]
			if key not in final_dictionary: # Must create the list
				new_value = [None, 0.0, None, 0.0, None, 0.0, 1.0] # use floats to ensure no flooring during division

				# Trees
				new_value[0] = value[0]
				new_value[2] = value[2]
				new_value[4] = value[4]

				# Support Values (P)
				new_value[1] = 1.0 if value[1] > bootstrap_cutoff else 0.0
				new_value[3] = 1.0 if value[3] > bootstrap_cutoff else 0.0
				new_value[5] = 1.0 if value[5] > bootstrap_cutoff else 0.0

				final_dictionary[key] = new_value

			else:
				# Update the support values (P)
				final_dictionary[key][1] += 1 if value[1] >= bootstrap_cutoff else 0
				final_dictionary[key][3] += 1 if value[3] >= bootstrap_cutoff else 0
				final_dictionary[key][5] += 1 if value[5] >= bootstrap_cutoff else 0

				# Update IQ spot to # of genes that have these 4 taxa present
				final_dictionary[key][6] += 1

	# Normalize support values (P) and add IQ values to the 6th index
	for key in final_dictionary:
		value = final_dictionary[key]
		value[1] = value[1] / value[6]
		value[3] = value[3] / value[6]
		value[5] = value[5] / value[6]

		iq = 1
		for index in [1, 3, 5]:
			if value[index] > 0:
				iq += (value[index] * log(value[index], 3)) # P(ti) * log base 3 of P(ti)

		value[6] = iq

		final_dictionary[key] = value

	return final_dictionary


# TEST getAllQuartetsAsDictionary
# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg3_5")
# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")
# trees = readTreeFile("testTree.txt")

# dictonary_of_quartets_1 = getAllQuartetsAsDictionary(trees)

# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")
# dictonary_of_quartets_3_5 = getAllQuartetsAsDictionary(trees)


# final_dictionary = mergeQuartetDictionaries([dictonary_of_quartets_1, dictonary_of_quartets_3_5], 8)

# counter = 0
# for key in final_dictionary:
# 	print "{0} - {1}: {2}".format(counter, key, final_dictionary[key])
# 	counter += 1




	# TEST getAllQuartetsAsDictionary
	# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg3_5")
	# trees = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")
	# trees = readTreeFile("testTree.txt")
	# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);

	# start_time = timeit.default_timer()
	# dictonary_of_quartets = getAllQuartetsAsDictionary(trees)
	# end_time = timeit.default_timer()

	# pprint(dictonary_of_quartets)
	# counter = 0
	# for key in dictonary_of_quartets:
	# 	print "{0} - {1}: {2}".format(counter, key, dictonary_of_quartets[key])
	# 	counter += 1

	# print "Dictionary Creation Time: {0}".format(end_time - start_time)


	# TEST putQuartetsDictionary
	# t = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")[0]
	# dictonary_of_quartets = combinationsOfLeaves(t)
	# putQuartetsDictionary(t, dictonary_of_quartets)
	# putQuartetsDictionary(t, to_prune, dictonary_of_quartets)
	# pprint(dictonary_of_quartets)
	# print len(dictonary_of_quartets)


	# TEST combinationsOfLeaves
	# t = readTreeFile("../for_isaac2/for_issac/RAxML_bootstrap.orfg1")[0]
	# print len(t)
	# pprint(len(combinationsOfLeaves(t)))



# input (list of branches, dictionary)
# for each branch:
#   generate all possible quartets (make bipartition, take all combinations of 2 for all combinations of 2 on the other side)
#   for each quartet:
#     if the quartet exists in the dictionary:
#       if the most frequent topology matches that of the reference tree -> map a 1
#       else -> -1
#       multiply the 1/-1 by the IC value to give it a weight
#    sum these values for all the quartets and place that value on the branch

def buildTree(referenceTreeFile, quartet_dictionary):
	# referenceTree = dendropy.Tree.get(path=referenceTreeFile, schema="newick")
	referenceTree = dendropy.Tree.get(data=referenceTreeFile, schema="newick")
	# print referenceTree
	tn = referenceTree.taxon_namespace
	print tn.labels()[::-1]

	splits = getListOfSplits(tn, referenceTree)
	# print splits, len(splits)

	# for split_object in splits:
	split_object = splits[0]
	left_combinations = list(combinations(split_object['left'], 2))
	right_combinations = list(combinations(split_object['right'], 2))
	# print(split_object['left'])
	# print(split_object['right'])
	# print(left_combinations)
	# print(right_combinations)


# List of objects with a right and a left
def getListOfSplits(taxonNamespace, tree):
	splits = []

	for bipartition in tree.encode_bipartitions():
		leafset = bipartition.leafset_taxa(taxonNamespace)
		# print len(leafset)

		# or len(leafset) == len(taxonNamespace) + 1 or len(leafset) == 1
		# if len(leafset) == len(taxonNamespace):
		# 	continue

		# tempTree = tree.clone()
		# tempTree2 = tree.clone()
		#
		# tempTree.prune_taxa(bipartition.leafset_taxa(taxonNamespace))
		# tempTree2.retain_taxa(bipartition.leafset_taxa(taxonNamespace))
		#
		# # Must have at least 2 leaves to form a quartet
		# if (len(tempTree.leaf_nodes()) < 2 or len(tempTree2.leaf_nodes()) < 2):
		# 	continue

		# print bipartition.split_as_bitstring() + " => " + bipartition.split_as_newick_string(taxonNamespace)
		split_object = getTaxaFromBipartition(taxonNamespace, bipartition)
		if split_object is not None:
			splits.append(split_object)

	return splits


# returns an object with 'left' and 'right' lists of taxa
# 1->left, 0->right
def getTaxaFromBipartition(taxonNamespace, bipartition):
	right = []
	left = []

	bString = bipartition.split_as_bitstring()

	index = len(taxonNamespace.labels()) - 1 # start at end because LSB
	for label in taxonNamespace.labels():
		if index < 0:
			break
		if bString[index] is '1':
			left.append(taxonNamespace.get_taxon(label))
			# left.append(label)
		else:
			right.append(taxonNamespace.get_taxon(label))
			# right.append(label)
		index -= 1

	return_object = {
		'right': right,
		'left': left
	}

	if len(return_object['left']) < 2 or len(return_object['right']) < 2:
		return None

	return return_object






s = "((A,B),(C,D));((A,C),(B,D));"
t = "((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);"
# tree will be '((A,B),(C,D))'
# buildTree("./for_issac/complete/RAxML_bestTree.orfg1.last_2", {})
buildTree(t, {})
