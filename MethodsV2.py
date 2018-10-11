#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from pprint import pprint


# Input: list of files with trees (requires newick format)
# Output: Dendropy TreeList object with all the trees
def readTrees(filenames):
    sample_tree_list = TreeList()
    for f in filenames:
        sample_tree_list.read(file=open(f, 'r'), schema="newick")
    return sample_tree_list

# Test readTrees
sample_tree_list = readTrees(['testTree.txt'])
# print(sample_tree_list.frequency_of_bipartition(labels=['H', 'E', 'C']))
# for t in sample_tree_list:
#     print(t, type(t.taxon_namespace[0]))




# MAKE DICT FOR EACH TREE
# COMBINE AT THE END



# Input: Dendropy TreeList object
# Output: Quartet dictionary with all unique quartets from the tree_list
def makeQuartetDictionary(tree_list):
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

# Test makeQuartetDictionary
quartet_dictionary = makeQuartetDictionary(sample_tree_list)
# pprint(len(quartet_dictionary))


# Input: Dendropy Tree object, a quartet dictionary as created by makeQuartetDictionary()
# Output: The quartet dictionary with updated support vectors for that tree
def getTreeQuartetSupport(tree, quartet_dictionary):
    taxon_label_list = [(n.taxon.label) for n in tree.leaf_nodes()]
    frozenset_of_taxa = frozenset(taxon_label_list) # unique set of all taxa

    print(tree.as_ascii_plot())

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
            result0 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[1]])
            if (result0 > 0):
                quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
                continue

            # Check 2nd Topology
            result1 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[2]])
            if (result1 > 0):
                quartet_dictionary[quartet][3] = quartet_dictionary[quartet][3] + 1
                continue

            # Check 3rd Topology
            result2 = single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[3]])
            if (result2 > 0):
                quartet_dictionary[quartet][5] = quartet_dictionary[quartet][5] + 1
                continue

            # ERROR
            print("Error: Topology is not a match")
            print(quartet, quartet_dictionary[quartet])
            break


getTreeQuartetSupport(sample_tree_list[0], quartet_dictionary)
# getTreeQuartetSupport(sample_tree_list[1], quartet_dictionary)
# pprint(quartet_dictionary)












# A - trees with duo (freq * number trees)
# B - trees with quad
# A >= B

# 100 trees total
# 50 trees have the quad
# 70 have the split


    # for quartet in quartet_dictionary:
    #     # print((quartet[0], quartet[1]), (quartet[0], quartet[2]), (quartet[0], quartet[3]))
    #     support_vector = quartet_dictionary[quartet]
    #
    #
    #
    #
    #
    #
    #
    #
    #
    #     print(quartet, support_vector[0], support_vector[2], support_vector[4])
    #     print(tree_list.frequency_of_bipartition(labels=[quartet[0], quartet[1]]))
    #     print(tree_list.frequency_of_bipartition(labels=[quartet[0], quartet[1], quartet[2], quartet[3]]))
        # for tree in tree_list:
        #     support_vector = quartet_dictionary[quartet]
        #     print(quartet, support_vector[0], support_vector[2], support_vector[4])


        #     pruned_tree = tree.extract_tree_with_taxa_labels(quartet)
        #     print(pruned_tree, support_vector[0], support_vector[2], support_vector[4])
        #     print(pruned_tree.print_plot())
        #     print(support_vector[0].print_plot())
        #     break
        #     if pruned_tree != None and support_vector != None:
        #         # Check 1st Topology
        #         result0 = treecompare.symmetric_difference(pruned_tree, support_vector[0])
        #         if (result0 is 0):
        #             quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
        #             continue
        #
        #         # Check 2nd Topology
        #         result1 = treecompare.symmetric_difference(pruned_tree, support_vector[2])
        #         if (result1 is 0):
        #             quartet_dictionary[quartet][3] = quartet_dictionary[quartet][3] + 1
        #             continue
        #
        #         # Check 3rd Topology
        #         result2 = treecompare.symmetric_difference(pruned_tree, support_vector[4])
        #         if (result2 is 0):
        #             quartet_dictionary[quartet][5] = quartet_dictionary[quartet][5] + 1
        #             continue
        #
        #         # ERROR
        #         # print("Error: Topology is not a match")
        #         # print(quartet, quartet_dictionary[quartet])
        #         # break
        # break
    # return quartet_dictionary

# Test fillQuartetSupportVectors
# quartet_dictionary = fillQuartetSupportVectors(sample_tree_list, quartet_dictionary)
# print(quartet_dictionary)














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
