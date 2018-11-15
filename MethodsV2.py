#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import json
import argparse
import sys
import os
import time
import multiprocessing


# STATS
# reference tree has 95 taxa
# smallest gene tree has 15, takes ~1min alone
# orfg10_5 has 54, (from 7+ hours to 53min!!! and only 30 minutes on the super computer)
# ~3 million quartets in reference tree (about same in full quartet dictionary)
# ~12 million quartets embedded in branches


# Input: list of files with trees (requires newick format)
# Output: A list of Dendropy TreeList objects with all the bootstrap trees for a gene tree in each
def readTrees(filenames, namespace, quiet=False):
    if not quiet:
        print()
        print("Reading in files...")
        print()

    sample_tree_list = []
    for f in filenames:
        # temp = TreeList(taxon_namespace=namespace)
        temp = TreeList()
        try:
            temp.read(file=open(f, 'r'), schema="newick", preserve_underscores=True)
        except:
            print("Error with file '{}': please only use files with newick tree format".format(f))
            sys.exit()

        sample_tree_list.append(temp)
    return sample_tree_list


# Input: Dendropy TreeList object
# Output: Quartet dictionary with all unique quartets from the tree_list
def makeQuartetDictionary(tree_list):
    taxon_label_list = tree_list.taxon_namespace.labels()
    combinations_of_taxa = combinations(taxon_label_list, 4)

    dictonary_of_quartets = {}

    for tuple_of_leaves in combinations_of_taxa:
        # sorted_list_of_leaves = list(tuple_of_leaves)
        # sorted_list_of_leaves.sort()
        frozenset_of_leaves = frozenset(tuple_of_leaves)
        dictonary_of_quartets[frozenset_of_leaves] = [0, 0, 0]

    return dictonary_of_quartets

# Input: Dendropy Tree object, a quartet dictionary as created by makeQuartetDictionary()
# Output: A new quartet dictionary with updated support vectors for that tree
def getTreeQuartetSupport(tree, quartet_dictionary, timing):
    tree.is_rooted = False
    tree.encode_bipartitions()

    taxon_label_list = tree.taxon_namespace.labels()

    frozenset_of_taxa = frozenset(taxon_label_list)  # unique set of all taxa
    bipartition_encoding = set(b.split_bitmask for b in tree.bipartition_encoding)
    bitstring_encoding = []
    for b in tree.bipartition_encoding:
        if not b.is_trivial():
            bitstring_encoding.append(b.split_as_bitstring())

    # pdm = tree.phylogenetic_distance_matrix()
    # pprint(pdm.distances(is_weighted_edge_distances=False))
    # sys.exit(1)

    bipartition_dictionary = makeBipartitionDictionary(taxon_label_list, bitstring_encoding)

    extraction_needed = 0
    counter = 0
    for quartet in quartet_dictionary:

        # if counter == 41:
        # if getShortestPath(quartet, pdm, tree):
        #     sys.exit(1)

        counter += 1
        if quartet.issubset(frozenset_of_taxa):  # if the tree contains the quartet
            p = round((counter / len(quartet_dictionary)) * 100, 2)
            e = round(extraction_needed/len(quartet_dictionary) * 100, 2)
            if timing:
                sys.stdout.write("                         Tree support progress: %f%% \t Extractions Needed: %f%%   \r" % (p, e) )
                sys.stdout.flush()
            try:
                dict_index = quartetBipartitionSupportHelper(tree, quartet_dictionary, quartet, bipartition_encoding, taxon_label_list, bitstring_encoding, bipartition_dictionary)
            # if dict_index < 0:
            #     sys.exit(1)
            except:
                extraction_needed += 1
                dict_index = quartetExtractionSupportHelper(tree, quartet_dictionary, quartet)
    return quartet_dictionary

def quartetBipartitionSupportHelper(tree, quartet_dictionary, quartet, bipartition_encoding, taxon_label_list, bitstring_encoding, bipartition_dictionary):
    sorted_quartet = list(quartet)
    sorted_quartet.sort()

    result0 = ((tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[1]]) in bipartition_encoding) or
              (tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[2], sorted_quartet[3]]) in bipartition_encoding))
    if (result0):
        quartet_dictionary[quartet][0] = quartet_dictionary[quartet][0] + 1
        return 0

    # Check 2nd Topology
    result1 = ((tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[2]]) in bipartition_encoding) or
              (tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[1], sorted_quartet[3]]) in bipartition_encoding))
    if (result1):
        quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
        return 1

    # Check 3rd Topology
    result2 = ((tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[3]]) in bipartition_encoding) or
              (tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[1], sorted_quartet[2]]) in bipartition_encoding))
    if (result2):
        quartet_dictionary[quartet][2] = quartet_dictionary[quartet][2] + 1
        return 2

    # dict_index = manualBitmaskSearch(sorted_quartet, taxon_label_list, bitstring_encoding)
    dict_index = manualBitmaskSearchV2(sorted_quartet, bipartition_dictionary)
    # compare = quartetExtractionSupportHelper(tree, quartet_dictionary, quartet)
    # if (dict_index is not compare and dict_index is not dict_indexV2):
        # dict_index = -1
    if (dict_index >= 0):
        quartet_dictionary[quartet][dict_index] = quartet_dictionary[quartet][dict_index] + 1
        return dict_index


    # ERROR
    if dict_index < 0:
        print()
        print("-----ERROR OUTPUT BIPARTITION-----")
        print("Sorted quartet:", sorted_quartet)
        print([taxon_label_list.index(sorted_quartet[0]), taxon_label_list.index(sorted_quartet[1]), taxon_label_list.index(sorted_quartet[2]), taxon_label_list.index(sorted_quartet[3])])
        print([len(taxon_label_list) - taxon_label_list.index(sorted_quartet[0]) - 1, len(taxon_label_list) - taxon_label_list.index(sorted_quartet[1]) - 1, len(taxon_label_list) - taxon_label_list.index(sorted_quartet[2]) - 1, len(taxon_label_list) - taxon_label_list.index(sorted_quartet[3]) - 1])
        print()
        print('Newick Tree: ', tree.as_string('newick'))
        print()
        print(tree.as_ascii_plot())
        print()
        print(taxon_label_list)
        print(tree.taxon_namespace.labels())
        print([[b.split_as_bitstring(), b.leafset_as_newick_string(tree.taxon_namespace)] for b in tree.bipartition_encoding])
        print()
        temp = [b.split_as_bitstring() for b in tree.bipartition_encoding]
        temp.sort()
        temp = [[t.taxon_bitmask] for t in quartet]
        print('temp ', temp)
        print(bitstring_encoding)
        print()
        print(tree.extract_tree_with_taxa_labels(quartet).as_ascii_plot())
        return dict_index

        # print(set(b.leafset_as_newick_string(tree.taxon_namespace) for b in tree.bipartition_encoding))
        # print(set(b.leafset_taxa(tree.taxon_namespace) for b in tree.bipartition_encoding))
        # print(tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[2]]))
        # print(manualBitmaskSearch(sorted_quartet[2], sorted_quartet[0], [(n.taxon.label) for n in tree.leaf_nodes()], tree.bipartition_encoding))
        # print(tree.extract_tree_with_taxa_labels(quartet).as_ascii_plot())
        pass
    raise Exception('Error: Topology is not a match')

def manualBitmaskSearch(sorted_quartet, taxon_label_list, bitstring_encoding):
    # Start from back because encoded LSB
    taxa_zero  = len(taxon_label_list) - taxon_label_list.index(sorted_quartet[0]) - 1
    taxa_one   = len(taxon_label_list) - taxon_label_list.index(sorted_quartet[1]) - 1
    taxa_two   = len(taxon_label_list) - taxon_label_list.index(sorted_quartet[2]) - 1
    taxa_three = len(taxon_label_list) - taxon_label_list.index(sorted_quartet[3]) - 1
    for b in bitstring_encoding:

        if (b[taxa_zero] is b[taxa_one]) and (b[taxa_two] is b[taxa_three]) and (b[taxa_zero] is not b[taxa_two]):
            return 0
        if (b[taxa_zero] is b[taxa_two]) and (b[taxa_one] is b[taxa_three]) and (b[taxa_zero] is not b[taxa_one]):
            return 1
        if (b[taxa_zero] is b[taxa_three]) and (b[taxa_one] is b[taxa_two]) and (b[taxa_zero] is not b[taxa_one]):
            return 2
    return -1

def manualBitmaskSearchV2(sorted_quartet, bipartition_dictionary):

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[1]])
    second_duet = frozenset([sorted_quartet[2], sorted_quartet[3]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 0
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 0
    # for b in bipartition_dictionary[first_duet]['1']:
    #     if b in bipartition_dictionary[second_duet]['0']:
    #         return 0
    # for b in bipartition_dictionary[second_duet]['1']:
    #     if b in bipartition_dictionary[first_duet]['0']:
    #         return 0

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[2]])
    second_duet = frozenset([sorted_quartet[1], sorted_quartet[3]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 1
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 1
    # for b in bipartition_dictionary[first_duet]['1']:
    #     if b in bipartition_dictionary[second_duet]['0']:
    #         return 1
    # for b in bipartition_dictionary[second_duet]['1']:
    #     if b in bipartition_dictionary[first_duet]['0']:
    #         return 1

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[3]])
    second_duet = frozenset([sorted_quartet[1], sorted_quartet[2]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 2
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 2
    # for b in bipartition_dictionary[first_duet]['1']:
    #     if b in bipartition_dictionary[second_duet]['0']:
    #         return 2
    # for b in bipartition_dictionary[second_duet]['1']:
    #     if b in bipartition_dictionary[first_duet]['0']:
    #         return 2
    return -1

def makeBipartitionDictionary(taxon_label_list, bitstring_encoding):
    # taxa duets
    combinations_of_taxa = combinations(taxon_label_list, 2)

    bipartition_dictionary = {}

    for tuple_of_leaves in combinations_of_taxa:
        # Start from back because encoded LSB
        taxa_zero  = len(taxon_label_list) - taxon_label_list.index(tuple_of_leaves[0]) - 1
        taxa_one   = len(taxon_label_list) - taxon_label_list.index(tuple_of_leaves[1]) - 1

        ones  = []
        zeros = []

        for b in bitstring_encoding:
            if (b[taxa_zero] is '1' and b[taxa_one] is '1'):
                ones.append(b)
            elif (b[taxa_zero] is '0' and b[taxa_one] is '0'):
                zeros.append(b)

        # frozenset_of_leaves = frozenset(tuple_of_leaves)
        bipartition_dictionary[frozenset(tuple_of_leaves)] = {
            '1': set(ones),
            '0': set(zeros)
        }
    return bipartition_dictionary


# Doesn't work
def getShortestPath(quartet, pdm, tree):
    sorted_quartet = list(quartet)
    sorted_quartet.sort()

    one = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[0]), tree.taxon_namespace.get_taxon(sorted_quartet[1]));
    two = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[2]), tree.taxon_namespace.get_taxon(sorted_quartet[3]));
    first = min(one, two);
    # print(first)

    three = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[0]), tree.taxon_namespace.get_taxon(sorted_quartet[2]));
    four = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[1]), tree.taxon_namespace.get_taxon(sorted_quartet[3]));
    second = min(three, four)
    # print(second)


    five = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[0]), tree.taxon_namespace.get_taxon(sorted_quartet[3]));
    six = pdm.path_edge_count(tree.taxon_namespace.get_taxon(sorted_quartet[1]), tree.taxon_namespace.get_taxon(sorted_quartet[2]));
    third = min(five, six)
    # print(third)

    total = min(first, second, third)
    if (total is first and total is second) or (total is first and total is third) or (total is second and total is third):
        print('quartet', sorted_quartet)
        print(tree.as_ascii_plot())
        print()
        print(one)
        print(two)
        print(three)
        print(four)
        print(five)
        print(six)
        print()
        print(tree.extract_tree_with_taxa_labels(quartet).as_ascii_plot())
        return True

    return False


def quartetExtractionSupportHelper(tree, quartet_dictionary, quartet):
    sorted_quartet = list(quartet)
    sorted_quartet.sort()

    # single_tree_list = TreeList()
    # single_tree_list.append(tree.extract_tree_with_taxa_labels(quartet))
    extracted_tree = tree.extract_tree_with_taxa_labels(quartet)
    extracted_tree.encode_bipartitions()
    bipartition_encoding = set(b.split_bitmask for b in extracted_tree.bipartition_encoding)

    # Check 1st Topology
    result0 = ((extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[1]]) in bipartition_encoding) or
              (extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[2], sorted_quartet[3]]) in bipartition_encoding))
    if (result0):
        quartet_dictionary[quartet][0] = quartet_dictionary[quartet][0] + 1
        return 0

    # Check 2nd Topology
    result1 = ((extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[2]]) in bipartition_encoding) or
              (extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[1], sorted_quartet[3]]) in bipartition_encoding))
    if (result1):
        quartet_dictionary[quartet][1] = quartet_dictionary[quartet][1] + 1
        return 1

    # Check 3rd Topology
    result2 = ((extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[0], sorted_quartet[3]]) in bipartition_encoding) or
              (extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[1], sorted_quartet[2]]) in bipartition_encoding))
    if (result2):
        quartet_dictionary[quartet][2] = quartet_dictionary[quartet][2] + 1
        return 2

    # ERROR
    # print("-----ERROR OUTPUT EXTRACTION-----")
    # print("Sorted quartet:", sorted_quartet)
    # print()
    # print(extracted_tree.as_string('newick'))
    # print()
    # print(extracted_tree.as_ascii_plot())
    # print()
    # print(bipartition_encoding)
    # print()
    # print(extracted_tree.taxon_namespace.taxa_bitmask(labels=[sorted_quartet[1], sorted_quartet[3]]))
    # print()
    raise Exception('Error: Topology is not a match')


# Input: An array of Dendropy TreeLists, a bootstrap cutoff value defaulting to 80
# Output: The full quartet dictionary for all gene trees containing P(t) and IC values
def buildFullSupport(gene_tree_list, bootstrap_cutoff_value=80, verbose=False, quiet=False, timing=False):
    if not quiet:
        print("Combining gene tree data into one dictionary...")
        print()
    full_quartet_dictionary = {}

    start = time.perf_counter()
    if timing:
        print('START: ', start)

    quartet_dictionary_queue = multiprocessing.Queue()
    processes = []
    quartet_dictionary_list = []

    for bootstrap_tree_list in gene_tree_list:
        quartet_dictionary_list.append(buildFullSupportParallelHelper(bootstrap_tree_list, start, verbose, quartet_dictionary_queue, timing))
    #     t = multiprocessing.Process(target=buildFullSupportParallelHelper, args=(bootstrap_tree_list, start, verbose, quartet_dictionary_queue))
    #     processes.append(t)
    #     t.start()
    # print("processes", processes)
    #
    # for one_process in processes:
    #     one_process.join()
    #
    # while not quartet_dictionary_queue.empty():
    #     quartet_dictionary_list.append(quartet_dictionary_queue.get())

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

def buildFullSupportParallelHelper(bootstrap_tree_list, start, verbose, parallel_queue, timing):
    quartet_dictionary = makeQuartetDictionary(bootstrap_tree_list)
    if timing:
        print('[PID: %d] MADE QUARTET DICTIONARY: ' % os.getpid(), (time.perf_counter() - start), '\t\t\tquartet_dictionary SIZE: ', len(quartet_dictionary))

    counter = 0
    for tree in bootstrap_tree_list:
        counter += 1
        getTreeQuartetSupport(tree, quartet_dictionary, timing)
        if timing:
            print('[%d/%d]' % (counter, len(bootstrap_tree_list)))
    if timing:
        print('[PID: %d] GOT FULL SUPPORT: ' % os.getpid(), (time.perf_counter() - start))

    if verbose:
        print("Full quartet dictionary:")
        [print(quartet, quartet_dictionary[quartet]) for quartet in quartet_dictionary]
        print()
    # parallel_queue.put(quartet_dictionary)
    return quartet_dictionary


def buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree="output_tree.tre", quiet=False, timing=False):
    reference_tree = Tree.get(path=referenceTreeFile, schema="newick", preserve_underscores=True)
    reference_tree.is_rooted = False
    reference_tree.encode_bipartitions()
    bipartition_encoding = set(b.split_bitmask for b in reference_tree.bipartition_encoding)

    bitstring_encoding = []
    for b in tree.bipartition_encoding:
        if not b.is_trivial():
            bitstring_encoding.append(b.split_as_bitstring())

    taxon_label_list = [(n.taxon.label) for n in reference_tree.leaf_nodes()]



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
        p = str(round((counter / len(splits)) * 100, 2)).rstrip('0')

        if 'left' not in split_object:
            continue


        left_combinations = list(combinations(split_object['left'], 2))
        right_combinations = list(combinations(split_object['right'], 2))

        total_possibilities = len(left_combinations) * len(right_combinations)

        total_support_value = 0

        inner_counter = 0
        for left_combination in left_combinations:
            for right_combination in right_combinations:
                inner_counter += 1
                inner_p = str(round((inner_counter / total_possibilities) * 100, 2)).rstrip('0')

                if timing:
                    sys.stdout.write("Build Label Tree Progress [%s/%s] : %s%%   %s%%\r" % (str(counter), str(len(splits)), p, inner_p))
                    sys.stdout.flush()

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


########## ARGPARSE
def runProgram(referenceTreeFile, sampleTreeList, bootstrap_cutoff_value=80, output_tree="output_tree.tre", verbose=False, quiet=False, timing=False):
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

    reference_tree_namespace = reference_tree.taxon_namespace

    sample_tree_list = readTrees(sampleTreeList, reference_tree_namespace, quiet)

    # Check if gene tree taxon namespace matches reference tree
    for s in sample_tree_list:
        if not reference_tree_namespace.has_taxa_labels(s.taxon_namespace.labels()):
            print('Error: reference tree is of a different taxon namespace as the sample trees')
            return

    full_quartet_dictionary = buildFullSupport(sample_tree_list, bootstrap_cutoff_value, verbose, quiet, timing)
    if verbose:
        print("Full quartet dictionary with support values")
        [print(quartet, full_quartet_dictionary[quartet])
               for quartet in full_quartet_dictionary]
        print()
    buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree, quiet, timing)

# ./MethodsV2.py run_files/RAxML_bestTree.rcGTA_cat run_files/RAxML_bootstrap.orfg1.last_2.subSample run_files/RAxML_bootstrap.orfg3_5.last_2.subSample -v -c 8
# ./MethodsV2.py test_trees/reference_tree.txt test_trees/trifurcations.txt -v -c 8
# ./MethodsV2.py test_trees/reference_tree.txt test_trees/highest_support.txt test_trees/highest_support.txt -v -c 8
# ./MethodsV2.py run_files/RAxML_bestTree.rcGTA_cat run_files/RAxML_bootstrap.orfg1.last_2 -v -c 8
# ./MethodsV2.py run_files/RAxML_bestTree.rcGTA_cat run_files/RAxML_bootstrap.orfg1.last_2 run_files/RAxML_bootstrap.orfg10.last_2 run_files/RAxML_bootstrap.orfg10_5.last_3 run_files/RAxML_bootstrap.orfg11.last_2 run_files/RAxML_bootstrap.orfg12.last_2 run_files/RAxML_bootstrap.orfg13.last_2 run_files/RAxML_bootstrap.orfg14.last_2 run_files/RAxML_bootstrap.orfg15.last_2 run_files/RAxML_bootstrap.orfg2.last_2 run_files/RAxML_bootstrap.orfg3.last_2 run_files/RAxML_bootstrap.orfg3_5.last_2 run_files/RAxML_bootstrap.orfg4.last_2 run_files/RAxML_bootstrap.orfg5.last_2 run_files/RAxML_bootstrap.orfg6.last_2 run_files/RAxML_bootstrap.orfg7.last_2 run_files/RAxML_bootstrap.orfg8.last_2 run_files/RAxML_bootstrap.orfg9.last_2 -v > run_output.txt
# ./MethodsV2.py run_files/RAxML_bestTree.rcGTA_cat run_files/RAxML_bootstrap.orfg1.last_2.subSample run_files/RAxML_bootstrap.orfg10.last_2.subSample run_files/RAxML_bootstrap.orfg10_5.last_3.subSample run_files/RAxML_bootstrap.orfg11.last_2.subSample run_files/RAxML_bootstrap.orfg12.last_2.subSample run_files/RAxML_bootstrap.orfg13.last_2.subSample run_files/RAxML_bootstrap.orfg14.last_2.subSample run_files/RAxML_bootstrap.orfg15.last_2.subSample run_files/RAxML_bootstrap.orfg2.last_2.subSample run_files/RAxML_bootstrap.orfg3.last_2.subSample run_files/RAxML_bootstrap.orfg3_5.last_2.subSample run_files/RAxML_bootstrap.orfg4.last_2.subSample run_files/RAxML_bootstrap.orfg5.last_2.subSample run_files/RAxML_bootstrap.orfg6.last_2.subSample run_files/RAxML_bootstrap.orfg7.last_2.subSample run_files/RAxML_bootstrap.orfg8.last_2.subSample run_files/RAxML_bootstrap.orfg9.last_2.subSample -v -c 8 > run_output.txt
parser = argparse.ArgumentParser()
parser.add_argument("reference_tree_file", metavar='<Reference Tree File>', help="The path of the reference tree file")
parser.add_argument('bootstrap_gene_tree_files', metavar='<Bootstrap Tree Files>', nargs='+',
                    help='The gene tree file paths containing bootstrap trees')
parser.add_argument("-v", "--verbose", action="store_true", default=False)
parser.add_argument("-q", "--quiet", action="store_true", default=False)
parser.add_argument("-t", "--timing", action="store_true", default=False, help='Setting timing turns verbose off.')
parser.add_argument("-c", "--cutoff", default=80, type=int,
                    help="Bootstrap Cutoff Value")
parser.add_argument("-o", "--output_file", default="output_tree.tre",
                    help="Output file for resulting tree with support. The default is 'output_tree.tre'")
args = parser.parse_args()

runProgram(args.reference_tree_file, args.bootstrap_gene_tree_files,
           bootstrap_cutoff_value=args.cutoff, output_tree=args.output_file, verbose=args.verbose, quiet=args.quiet, timing=args.timing)
