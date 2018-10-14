#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import json


# Input: list of files with trees (requires newick format)
# Output: A list of Dendropy TreeList objects with all the bootstrap trees for a gene tree in each
def readTrees(filenames):
    print("Reading in files...")
    sample_tree_list = []
    for f in filenames:
        temp = TreeList()
        temp.read(file=open(f, 'r'), schema="newick", preserve_underscores=True)
        sample_tree_list.append(temp)
    return sample_tree_list

# TEST readTrees
# sample_tree_list = readTrees(['testTree.txt'])
sample_tree_list = readTrees(['for_issac/complete/RAxML_bootstrap.orfg1.last_2'])
# 'for_issac/complete/RAxML_bootstrap.orfg3_5.last_2'

# Input: Dendropy TreeList object
# Output: Quartet dictionary with all unique quartets from the tree_list
def makeQuartetDictionary(tree_list):
    print("Making the quartet dictionary...")
    taxon_label_list = tree_list.taxon_namespace.labels()
    combinations_of_taxa = combinations(taxon_label_list, 4)

    dictonary_of_quartets = {}

    for tuple_of_leaves in combinations_of_taxa:
        sorted_list_of_leaves = list(tuple_of_leaves)
        sorted_list_of_leaves.sort()
        frozenset_of_leaves = frozenset(sorted_list_of_leaves)
        topology1 = "(({0},{1}),({2},{3}));".format(sorted_list_of_leaves[0], sorted_list_of_leaves[1], sorted_list_of_leaves[2], sorted_list_of_leaves[3])
        topology2 = "(({0},{1}),({2},{3}));".format(sorted_list_of_leaves[0], sorted_list_of_leaves[2], sorted_list_of_leaves[1], sorted_list_of_leaves[3])
        topology3 = "(({0},{1}),({2},{3}));".format(sorted_list_of_leaves[0], sorted_list_of_leaves[3], sorted_list_of_leaves[1], sorted_list_of_leaves[2])
        dictonary_of_quartets[frozenset_of_leaves] = [Tree.get(data=topology1, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0,
                                                  Tree.get(data=topology2, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0,
                                                  Tree.get(data=topology3, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0]

    return dictonary_of_quartets

# TEST makeQuartetDictionary
# quartet_dictionary = makeQuartetDictionary(sample_tree_list)
# pprint(len(quartet_dictionary))


# Input: Dendropy Tree object, a quartet dictionary as created by makeQuartetDictionary()
# Output: A new quartet dictionary with updated support vectors for that tree
def getTreeQuartetSupport(tree, quartet_dictionary):
    taxon_label_list = [(n.taxon.label) for n in tree.leaf_nodes()]
    frozenset_of_taxa = frozenset(taxon_label_list) # unique set of all taxa

    for quartet in quartet_dictionary:
        if quartet.issubset(frozenset_of_taxa): # if the tree contains the quartet
            sorted_quartet = list(quartet)
            sorted_quartet.sort()

            single_tree_list = TreeList()
            single_tree_list.append(tree.extract_tree_with_taxa_labels(quartet))

            # print(sorted_quartet, single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[1]]),
            #                       single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[2]]),
            #                       single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[3]]))

            # Check 1st Topology
            result0 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[1]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[2], sorted_quartet[3]])
            # print(result0)
            if (result0 > 0):
                quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
                continue

            # Check 2nd Topology
            result1 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[2]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[1], sorted_quartet[3]])
            if (result1 > 0):
                quartet_dictionary[quartet][3] = quartet_dictionary[quartet][3] + 1
                continue

            # Check 3rd Topology
            result2 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[3]])  + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[1], sorted_quartet[2]])
            if (result2 > 0):
                quartet_dictionary[quartet][5] = quartet_dictionary[quartet][5] + 1
                continue

            # ERROR
            # print(single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[1]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[2], sorted_quartet[3]]))
            print("Error: Topology is not a match")
            print()
            print("Sorted quartet:", sorted_quartet)
            print("\n")
            print(single_tree_list[0].as_ascii_plot())
            quit()

    return quartet_dictionary

# TEST getTreeQuartetSupport
# getTreeQuartetSupport(sample_tree_list[0], quartet_dictionary)
# getTreeQuartetSupport(sample_tree_list[1], quartet_dictionary)
# pprint(quartet_dictionary)


# Input: An array of Dendropy TreeLists, a bootstrap cutoff value defaulting to 80
# Output: The full quartet dictionary for all gene trees containing P(t) and IC values
def buildFullSupport(gene_tree_list, bootstrap_cutoff_value=80):
    print("Combining gene tree data into one dictionary...")
    full_quartet_dictionary = {}

    for bootstrap_tree_list in gene_tree_list:
        quartet_dictionary = makeQuartetDictionary(bootstrap_tree_list)
        for tree in bootstrap_tree_list:
            getTreeQuartetSupport(tree, quartet_dictionary)
        # [print(quartet, quartet_dictionary[quartet]) for quartet in quartet_dictionary]

        # Find support > 80 and add a 1 in the full dictionary, otherwise add a 0
        # Use index 6 of the list to record how many times the quartet is seen
        for quartet in quartet_dictionary:
            if quartet not in full_quartet_dictionary:
                full_quartet_dictionary[quartet] = [quartet_dictionary[quartet][0], 0.0, quartet_dictionary[quartet][2], 0.0, quartet_dictionary[quartet][4], 0.0, 1.0]
            else:
                full_quartet_dictionary[quartet][6] += 1.0
            full_quartet_dictionary[quartet][1] += (1.0 if quartet_dictionary[quartet][1] > bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][3] += (1.0 if quartet_dictionary[quartet][3] > bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][5] += (1.0 if quartet_dictionary[quartet][5] > bootstrap_cutoff_value else 0.0)

    # full_quartet_dictionary now has all the support vectors, s(t)
    # Next we must normalize the support vectors: p(t) = s(ti)/(s(t1) + s(t2) + s(t3))
    # Also generate the IQ value
    for quartet in full_quartet_dictionary:
        # Normalize the support values
        full_quartet_dictionary[quartet][1] = full_quartet_dictionary[quartet][1] / full_quartet_dictionary[quartet][6]
        full_quartet_dictionary[quartet][3] = full_quartet_dictionary[quartet][3] / full_quartet_dictionary[quartet][6]
        full_quartet_dictionary[quartet][5] = full_quartet_dictionary[quartet][5] / full_quartet_dictionary[quartet][6]

        iq = 1.0
        for index in [1, 3, 5]:
            if full_quartet_dictionary[quartet][index] > 0:
                iq += (full_quartet_dictionary[quartet][index] * log(full_quartet_dictionary[quartet][index], 3)) # P(ti) * log base 3 of P(ti)

        full_quartet_dictionary[quartet][6] = iq

    return full_quartet_dictionary

# TEST buildFullSupport
full_quartet_dictionary = buildFullSupport(sample_tree_list)
[print(quartet, full_quartet_dictionary[quartet]) for quartet in full_quartet_dictionary]










# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((D,((C,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((C,((D,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
# ((A,B),((D,((C,(E,(F,G))),H)),(((I,(J,K)),(L,M)),N)),O);
