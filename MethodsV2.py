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


# Input: Dendropy TreeList object
# Output: Quartet dictionary with all unique quartets from the tree_list
def makeQuartetDictionary(tree_list):
    taxon_list = tree_list.taxon_namespace.get_taxa
    combinations_of_taxa = combinations(taxon_list)
    
