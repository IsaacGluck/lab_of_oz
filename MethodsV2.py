#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import json
import argparse


# Input: list of files with trees (requires newick format)
# Output: A list of Dendropy TreeList objects with all the bootstrap trees for a gene tree in each
def readTrees(filenames, quiet=False):
    if not quiet:
        print()
        print("Reading in files...")
        print()
    sample_tree_list = []
    for f in filenames:
        temp = TreeList()
        temp.read(file=open(f, 'r'), schema="newick",
                  preserve_underscores=True)
        sample_tree_list.append(temp)
    return sample_tree_list


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
        topology1 = "(({0},{1}),({2},{3}));".format(
            sorted_list_of_leaves[0], sorted_list_of_leaves[1], sorted_list_of_leaves[2], sorted_list_of_leaves[3])
        topology2 = "(({0},{1}),({2},{3}));".format(
            sorted_list_of_leaves[0], sorted_list_of_leaves[2], sorted_list_of_leaves[1], sorted_list_of_leaves[3])
        topology3 = "(({0},{1}),({2},{3}));".format(
            sorted_list_of_leaves[0], sorted_list_of_leaves[3], sorted_list_of_leaves[1], sorted_list_of_leaves[2])
        dictonary_of_quartets[frozenset_of_leaves] = [Tree.get(data=topology1, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0,
                                                  Tree.get(
                                                      data=topology2, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0,
                                                  Tree.get(data=topology3, taxon_namespace=tree_list.taxon_namespace, schema="newick"), 0]

    return dictonary_of_quartets


# Input: Dendropy Tree object, a quartet dictionary as created by makeQuartetDictionary()
# Output: A new quartet dictionary with updated support vectors for that tree
def getTreeQuartetSupport(tree, quartet_dictionary):
    taxon_label_list = [(n.taxon.label) for n in tree.leaf_nodes()]
    frozenset_of_taxa = frozenset(taxon_label_list)  # unique set of all taxa

    for quartet in quartet_dictionary:
        if quartet.issubset(frozenset_of_taxa):  # if the tree contains the quartet
            sorted_quartet = list(quartet)
            sorted_quartet.sort()

            single_tree_list = TreeList()
            single_tree_list.append(
                tree.extract_tree_with_taxa_labels(quartet))

            # print(sorted_quartet, single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[x], sorted_quartet[1]]),
            #                       single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[2]]),
            #                       single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[0], sorted_quartet[3]]))

            # Check 1st Topology
            result0 = single_tree_list.frequency_of_bipartition(
                labels=[sorted_quartet[0], sorted_quartet[1]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[2], sorted_quartet[3]])
            # print(result0)
            if (result0 > 0):
                quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
                continue

            # Check 2nd Topology
            result1 = single_tree_list.frequency_of_bipartition(
                labels=[sorted_quartet[0], sorted_quartet[2]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[1], sorted_quartet[3]])
            if (result1 > 0):
                quartet_dictionary[quartet][3] = quartet_dictionary[quartet][3] + 1
                continue

            # Check 3rd Topology
            result2 = single_tree_list.frequency_of_bipartition(
                labels=[sorted_quartet[0], sorted_quartet[3]]) + single_tree_list.frequency_of_bipartition(labels=[sorted_quartet[1], sorted_quartet[2]])
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


# Input: An array of Dendropy TreeLists, a bootstrap cutoff value defaulting to 80
# Output: The full quartet dictionary for all gene trees containing P(t) and IC values
def buildFullSupport(gene_tree_list, bootstrap_cutoff_value=80, verbose=False, quiet=False):
    if not quiet:
        print("Combining gene tree data into one dictionary...")
        print()
    full_quartet_dictionary = {}

    for bootstrap_tree_list in gene_tree_list:
        quartet_dictionary = makeQuartetDictionary(bootstrap_tree_list)
        for tree in bootstrap_tree_list:
            getTreeQuartetSupport(tree, quartet_dictionary)

        if verbose:
            print("Full quartet dictionary:")
            [print(quartet, quartet_dictionary[quartet]) for quartet in quartet_dictionary]
            print()

        # Find support > 80 and add a 1 in the full dictionary, otherwise add a 0
        # Use index 6 of the list to record how many times the quartet is seen
        for quartet in quartet_dictionary:
            if quartet not in full_quartet_dictionary:
                full_quartet_dictionary[quartet] = [quartet_dictionary[quartet][0], 0.0,
                    quartet_dictionary[quartet][2], 0.0, quartet_dictionary[quartet][4], 0.0, 1.0]
            else:
                full_quartet_dictionary[quartet][6] += 1.0
            full_quartet_dictionary[quartet][1] += (
                1.0 if quartet_dictionary[quartet][1] >= bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][3] += (
                1.0 if quartet_dictionary[quartet][3] >= bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][5] += (
                1.0 if quartet_dictionary[quartet][5] >= bootstrap_cutoff_value else 0.0)

    # Probably don't need
    # print("Full quartet dictionary with support values MIDDLE")
    # [print(quartet, full_quartet_dictionary[quartet])
    #        for quartet in full_quartet_dictionary]
    # print()

    # full_quartet_dictionary now has all the support vectors, s(t)
    # Next we must normalize the support vectors: p(t) = s(ti)/(s(t1) + s(t2) + s(t3))
    # Also generate the IQ value
    for quartet in full_quartet_dictionary:
        # Normalize the support values
        full_quartet_dictionary[quartet][1] = full_quartet_dictionary[quartet][1] / \
            full_quartet_dictionary[quartet][6]
        full_quartet_dictionary[quartet][3] = full_quartet_dictionary[quartet][3] / \
            full_quartet_dictionary[quartet][6]
        full_quartet_dictionary[quartet][5] = full_quartet_dictionary[quartet][5] / \
            full_quartet_dictionary[quartet][6]

        iq = 1.0
        for index in [1, 3, 5]:
            if full_quartet_dictionary[quartet][index] > 0:
                # P(ti) * log base 3 of P(ti)
                iq += (full_quartet_dictionary[quartet][index] *
                       log(full_quartet_dictionary[quartet][index], 3))

        full_quartet_dictionary[quartet][6] = iq

    return full_quartet_dictionary


def buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree="output_tree.tre", quiet=False):
    reference_tree = Tree.get(path=referenceTreeFile, schema="newick")
    reference_tree.is_rooted = True
    reference_tree_list = TreeList()
    reference_tree_list.append(reference_tree)

    tn = reference_tree.taxon_namespace

    splits = getListOfSplits(tn, reference_tree)

    for split_object in splits:
        if 'left' not in split_object:
            continue
        left_combinations = list(combinations(split_object['left'], 2))
        right_combinations = list(combinations(split_object['right'], 2))

        total_possibilities = len(left_combinations) * len(right_combinations)

        total_support_value = 0

        for left_combination in left_combinations:
            for right_combination in right_combinations:
                combined_taxa = list(left_combination + right_combination)
                combined_taxa_labels = [taxa.label for taxa in combined_taxa]
                combined_taxa_labels.sort()
                quartet_dictionary_key = frozenset(combined_taxa_labels)
                if quartet_dictionary_key in full_quartet_dictionary:
                    extracted_tree = reference_tree.extract_tree_with_taxa(combined_taxa)
                    quartet_dictionary_value = full_quartet_dictionary[quartet_dictionary_key]

                    # indices of tree structures in dictionary
                    support_value = -1
                    result0 = reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[0], combined_taxa_labels[1]]) + reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[2], combined_taxa_labels[3]])
                    result1 = reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[0], combined_taxa_labels[2]]) + reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[1], combined_taxa_labels[3]])
                    result2 = reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[0], combined_taxa_labels[3]]) + reference_tree_list.frequency_of_bipartition(labels=[combined_taxa_labels[1], combined_taxa_labels[2]])
                    results = [result0, result1, result2]
                    for i in range(3):
                        if results[i] > 0:
                            max_val = max(quartet_dictionary_value[1], quartet_dictionary_value[3], quartet_dictionary_value[5])
                            if max_val == quartet_dictionary_value[(i * 2) + 1] and max_val != 0:
                                support_value = 1
                    support_value *= quartet_dictionary_value[6]
                    total_support_value += support_value
        split_object['edge'].head_node.label = total_support_value / total_possibilities
    if not quiet:
        print(reference_tree.as_string(schema="newick"))
    reference_tree.write(path=output_tree, schema="newick")


def getListOfSplits(taxonNamespace, tree):
    splits = []

    tree.encode_bipartitions()

    for node in tree:
        split_object = getTaxaFromBipartition(
            taxonNamespace, node.edge.bipartition)
        # if split_object is not None:
        split_object['edge'] = node.edge
        splits.append(split_object)
        # print(split_object)

    return splits

# returns an object with 'left' and 'right' lists of taxa
# 1->left, 0->right
def getTaxaFromBipartition(taxonNamespace, bipartition):
    right, left = [], []

    bString = bipartition.split_as_bitstring()

    index = len(taxonNamespace.labels()) - 1  # start at end because LSB
    for label in taxonNamespace.labels():
        if index < 0:
            break
        if bString[index] is '1':
            left.append(taxonNamespace.get_taxon(label))
        else:
            right.append(taxonNamespace.get_taxon(label))
        index -= 1

    return_object = {
        'right': right,
        'left': left
    }

    if len(return_object['left']) < 2 or len(return_object['right']) < 2:
        return {}

    return return_object


def runProgram(referenceTreeFile, sampleTreeList, bootstrap_cutoff_value=80, output_tree="output_tree.tre", verbose=False, quiet=False):
    if verbose:
        print("Reference Tree: ", referenceTreeFile)
        print("Sample Tree List: ", sampleTreeList)
        print("Bootstrap Cutoff Value: ", bootstrap_cutoff_value)
        print("Output Tree File: ", output_tree)
    sample_tree_list = readTrees(sampleTreeList, quiet)
    full_quartet_dictionary = buildFullSupport(sample_tree_list, bootstrap_cutoff_value, verbose, quiet)
    if verbose:
        print("Full quartet dictionary with support values")
        [print(quartet, full_quartet_dictionary[quartet])
               for quartet in full_quartet_dictionary]
        print()
    buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree, quiet)


parser = argparse.ArgumentParser()
parser.add_argument("reference_tree_file", metavar='Reference Tree Files', help="The path of the reference tree file")
parser.add_argument('bootstrap_gene_tree_files', metavar='Bootstrap Tree Files', nargs='+',
                    help='The gene tree file paths containing bootstrap trees')
parser.add_argument("-v", "--verbose", action="store_true", default=False)
parser.add_argument("-q", "--quiet", action="store_true", default=False)
parser.add_argument("-c", "--cutoff", default=80, type=int,
                    help="Bootstrap Cutoff Value")
parser.add_argument("-o", "--output_file", default="output_tree.tre",
                    help="Output file for resulting tree with support")
args = parser.parse_args()
# runProgram(referenceTreeFile, sampleTreeList, bootstrap_cutoff_value=80, output_tree="output_tree.tre", verbose=False, quiet=False):
runProgram(args.reference_tree_file, args.bootstrap_gene_tree_files,
           bootstrap_cutoff_value=args.cutoff, output_tree=args.output_file, verbose=args.verbose, quiet=args.quiet)
