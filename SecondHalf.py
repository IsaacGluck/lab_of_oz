#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import pickle, sys, os, time


# Input: list of files with serialized quartet_dictionaries
# Output: A list of quartet_dictionary objects
def readPickledTrees(quartetDictionaryFileList, quiet):
    if not quiet:
        print()
        print("Reading in quartet dictionary files...")
        print()

    quartet_dictionary_list = []
    for filename in quartetDictionaryFileList:
        file_obj = open(filename, 'rb')
        qd = pickle.load(file_obj)
        quartet_dictionary_list.append(qd)
        file_obj.close()

    return quartet_dictionary_list


# Input: An array of Dendropy TreeLists, a bootstrap cutoff value defaulting to 80
# Output: The full quartet dictionary for all gene trees containing P(t) and IC values
def buildFullSupport(quartet_dictionary_list, bootstrap_cutoff_value=80, verbose=False, quiet=False, timing=False):
    if not quiet:
        print("Combining gene tree data into one dictionary...")
        print()

    full_quartet_dictionary = {}

    start = time.perf_counter()
    if timing:
        print('START: ', start)

    for quartet_dictionary in quartet_dictionary_list:
        # Find support > 80 and add a 1 in the full dictionary, otherwise add a 0
        # Use index 6 of the list to record how many times the quartet is seen
        for quartet in quartet_dictionary:
            if quartet not in full_quartet_dictionary:
                full_quartet_dictionary[quartet] = [0.0, 0.0, 0.0, 1.0]
            else:
                full_quartet_dictionary[quartet][3] += 1.0
            full_quartet_dictionary[quartet][0] += (
                1.0 if quartet_dictionary[quartet][0] >= bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][1] += (
                1.0 if quartet_dictionary[quartet][1] >= bootstrap_cutoff_value else 0.0)
            full_quartet_dictionary[quartet][2] += (
                1.0 if quartet_dictionary[quartet][2] >= bootstrap_cutoff_value else 0.0)
        if timing:
            print('CONVERTED SUPPORT ABOVE CUTOFF TO 1: ', (time.perf_counter() - start), '\t\t\tfull_quartet_dictionary SIZE: ', len(full_quartet_dictionary))

    # full_quartet_dictionary now has all the support vectors, s(t)
    # Next we must normalize the support vectors: p(t) = s(ti)/(s(t1) + s(t2) + s(t3))
    # Also generate the IQ value
    for quartet in full_quartet_dictionary:
        # Normalize the support values
        full_quartet_dictionary[quartet][0] = full_quartet_dictionary[quartet][0] / \
            full_quartet_dictionary[quartet][3]
        full_quartet_dictionary[quartet][1] = full_quartet_dictionary[quartet][1] / \
            full_quartet_dictionary[quartet][3]
        full_quartet_dictionary[quartet][2] = full_quartet_dictionary[quartet][2] / \
            full_quartet_dictionary[quartet][3]

        iq = 1.0
        for index in range(3):
            if full_quartet_dictionary[quartet][index] > 0:
                # P(ti) * log base 3 of P(ti)
                iq += (full_quartet_dictionary[quartet][index] *
                       log(full_quartet_dictionary[quartet][index], 3))

        full_quartet_dictionary[quartet][3] = iq

    if timing:
        print('IC VALUES COMPUTED: ', (time.perf_counter() - start))
        print()

    return full_quartet_dictionary


def buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree="output_tree.tre", quiet=False, timing=False):
    reference_tree = Tree.get(path=referenceTreeFile, schema="newick", preserve_underscores=True)
    reference_tree.is_rooted = True
    reference_tree.encode_bipartitions()
    bipartition_encoding = set(b.split_bitmask for b in reference_tree.bipartition_encoding)
    taxon_label_list = [(n.taxon.label) for n in reference_tree.leaf_nodes()]
    # reference_tree_list = TreeList()
    # reference_tree_list.append(reference_tree)


    tn = reference_tree.taxon_namespace

    start = start = time.perf_counter()
    if timing:
        print('START BUILD LABEL TREE: ', start)
    splits = getListOfSplits(tn, reference_tree)
    if timing:
        print('GOT LIST OF SPLITS: ', (time.perf_counter() - start))




    counter = 0
    for split_object in splits:
        counter += 1
        p = round((counter / len(splits)) * 100, 2)
        if timing:
            sys.stdout.write("Build Label Tree Progress: %f%%   \r" % (p) )
            sys.stdout.flush()

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
                    quartet_dictionary_value = full_quartet_dictionary[quartet_dictionary_key]

                    # indices of tree structures in dictionary
                    support_value = -1
                    result0 = ((reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[0], combined_taxa_labels[1]]) in bipartition_encoding) or
                              (reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[2], combined_taxa_labels[3]]) in bipartition_encoding))
                    result1 = ((reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[0], combined_taxa_labels[2]]) in bipartition_encoding) or
                              (reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[1], combined_taxa_labels[3]]) in bipartition_encoding))
                    result2 = ((reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[0], combined_taxa_labels[3]]) in bipartition_encoding) or
                              (reference_tree.taxon_namespace.taxa_bitmask(labels=[combined_taxa_labels[1], combined_taxa_labels[2]]) in bipartition_encoding))
                    results = [result0, result1, result2]
                    for i in range(3):
                        if results[i]:
                            max_val = max(quartet_dictionary_value[0], quartet_dictionary_value[1], quartet_dictionary_value[2])
                            if max_val == quartet_dictionary_value[i] and max_val != 0:
                                support_value = 1
                    support_value *= quartet_dictionary_value[3]
                    total_support_value += support_value
        split_object['edge'].head_node.label = total_support_value / total_possibilities
    if not quiet:
        print(reference_tree.as_string(schema="newick", suppress_internal_node_labels=False))

    if timing:
        print('REFERENCE TREE BUILT: ', (time.perf_counter() - start))
        print()
    reference_tree.write(path=output_tree, schema="newick", suppress_internal_node_labels=False)

def getListOfSplits(taxonNamespace, tree):
    splits = []

    tree.encode_bipartitions()

    for node in tree:
        split_object = getTaxaFromBipartition(
            taxonNamespace, node.edge.bipartition)
        split_object['edge'] = node.edge
        splits.append(split_object)

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


def runProgram(referenceTreeFile, quartetDictionaryFileList, bootstrap_cutoff_value=80, output_tree="output_tree.tre", verbose=False, quiet=False, timing=False):
    if verbose:
        print("Reference Tree: ", referenceTreeFile)
        print("Sample Tree List: ", sampleTreeList)
        print("Bootstrap Cutoff Value: ", bootstrap_cutoff_value)
        print("Output Tree File: ", output_tree)

    if timing:
        verbose = False

    try:
        reference_tree = Tree.get(path=referenceTreeFile, schema="newick", preserve_underscores=True)
    except:
        print("Error with file '{}': please only use files with newick tree format".format(referenceTreeFile))
        sys.exit()

    quartet_dictionary_list = readPickledTrees(quartetDictionaryFileList, quiet)

    full_quartet_dictionary = buildFullSupport(quartet_dictionary_list, bootstrap_cutoff_value, verbose, quiet, timing)
    if verbose:
        print("Full quartet dictionary with support values")
        [print(quartet, full_quartet_dictionary[quartet])
               for quartet in full_quartet_dictionary]
        print()
    buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree, quiet, timing)

quartetDictionaryFileList = ['quartet_dictionaries/highest_support.txt.quartet_dictionary', 'quartet_dictionaries/low_support.txt.quartet_dictionary', 'quartet_dictionaries/medium_support.txt.quartet_dictionary']
# ./MethodsV2.py test_trees/reference_tree.txt test_trees/highest_support.txt test_trees/low_support.txt test_trees/medium_support.txt -t -c 8
runProgram('test_trees/reference_tree.txt', quartetDictionaryFileList, bootstrap_cutoff_value=8, timing=True)