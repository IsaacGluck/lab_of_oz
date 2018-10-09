#!/usr/bin/env python
import dendropy
from itertools import combinations
from dendropy.calculate import treecompare

# Input: list of files with trees (requires newick format)
# Output: Dendropy TreeList object with all the trees
def readTrees(filenames):
    sample_tree_list = TreeList()
    for f in filenames:
        sample_tree_list.read(file=open(f, 'r'), schema="newick")
    return sample_tree_list

# Test
sample_tree_list = readTrees(['testTree.txt'])
print(sample_tree_list)

# Input: Dendropy TreeList object
# Output: Quartet dictionary with all unique quartets from the tree_list
def makeQuartetDictionary(tree_list):
    taxon_list = tree_list.taxon_namespace.get_taxa
    combinations_of_taxa = combinations(taxon_list)

    dictonary_of_quartets = {}

    for tuple_of_leaves in combinations_of_taxa:
        topology1 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[1], tuple_of_leaf_names[2], tuple_of_leaf_names[3])
		topology2 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[2], tuple_of_leaf_names[1], tuple_of_leaf_names[3])
		topology3 = "(({0},{1}),{2},{3});".format(tuple_of_leaf_names[0], tuple_of_leaf_names[3], tuple_of_leaf_names[1], tuple_of_leaf_names[2])
        dictonary_of_quartets[tuple_of_leaves] = [Tree.get(data=topology1, schema="newick"), 0,
                                                  Tree.get(data=topology2, schema="newick"), 0,
                                                  Tree.get(data=topology3, schema="newick"), 0]

    return dictonary_of_quartets
