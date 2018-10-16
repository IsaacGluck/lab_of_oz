#!/usr/bin/env python
from MethodsV2 import *

# sumtrees.py --decimals=0 --percentages --output-tree-filepath=result.tre highest_support.txt

# TEST readTrees
# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/highest_support.txt'])
# sample_tree_list = readTrees(['test_trees/high_conflict.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/low_support.txt'])
# sample_tree_list = readTrees(['test_trees/medium_support.txt', 'test_trees/medium_support.txt'])

# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/low_support.txt'])
# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/medium_support.txt'])
# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/medium_support.txt'])

# sample_tree_list = readTrees(['for_issac/complete/RAxML_bootstrap.orfg1.last_2'])


# TEST makeQuartetDictionary
# quartet_dictionary = makeQuartetDictionary(sample_tree_list)
# pprint(len(quartet_dictionary))


# TEST getTreeQuartetSupport
# getTreeQuartetSupport(sample_tree_list[0], quartet_dictionary)
# getTreeQuartetSupport(sample_tree_list[1], quartet_dictionary)
# pprint(quartet_dictionary)


# TEST buildFullSupport
# full_quartet_dictionary = buildFullSupport(sample_tree_list, 8)
# print("Full quartet dictionary with support values")
# [print(quartet, full_quartet_dictionary[quartet])
#        for quartet in full_quartet_dictionary]
# print()

# buildLabeledTree("test_trees/reference_tree.txt", full_quartet_dictionary)
# buildLabeledTree("for_issac/complete/RAxML_bestTree.rcGTA_cat", full_quartet_dictionary)


# TEST buildLabeledTree
# buildLabeledTree("test_trees/reference_tree.txt", full_quartet_dictionary)

# TEST runProgram
# runProgram(referenceTreeFile, sampleTreeList, bootstrap_cutoff_value=80, output_tree="output_tree.tre", verbose=False, quiet=False):
runProgram("test_trees/reference_tree.txt", ['test_trees/highest_support.txt', 'test_trees/highest_support.txt'], 8, verbose=True, quiet=True)

# Full Test
# Reference Tree: for_issac/complete/RAxML_bestTree.rcGTA_cat
# all 17 bootstrap trees
