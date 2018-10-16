#!/usr/bin/env python
from MethodsV2 import *

# sumtrees.py --decimals=0 --percentages --output-tree-filepath=result.tre highest_support.txt

# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/highest_support.txt'])
# sample_tree_list = readTrees(['test_trees/high_conflict.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/low_support.txt'])
sample_tree_list = readTrees(['test_trees/medium_support.txt', 'test_trees/medium_support.txt'])


# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/low_support.txt'])
# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/medium_support.txt'])
# sample_tree_list = readTrees(['test_trees/highest_support.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/high_conflict.txt'])
# sample_tree_list = readTrees(['test_trees/low_support.txt', 'test_trees/medium_support.txt'])



full_quartet_dictionary = buildFullSupport(sample_tree_list, 8)
print("Full quartet dictionary with support values")
[print(quartet, full_quartet_dictionary[quartet])
       for quartet in full_quartet_dictionary]
print()

buildLabeledTree("test_trees/reference_tree.txt", full_quartet_dictionary)




# Full Test
# Reference Tree: for_issac/complete/RAxML_bestTree.rcGTA_cat
# all 17 bootstrap trees
