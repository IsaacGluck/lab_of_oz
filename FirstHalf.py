#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import pickle, sys, os, time


# Input: Filename with trees (requires newick format)
# Output: A Dendropy TreeList object with all the bootstrap trees for the gene tree
def readTree(filename, quiet=False):
    if not quiet:
        print()
        print("Reading in files...")
        print()

    temp = TreeList()
    try:
        temp.read(file=open(filename, 'r'), schema="newick", preserve_underscores=True)
    except:
        print("Error with file '{}': please only use files with newick tree format".format(f))
        sys.exit()

    return temp


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

    bipartition_dictionary = makeBipartitionDictionary(taxon_label_list, bitstring_encoding)

    extraction_needed = 0
    counter = 0
    for quartet in quartet_dictionary:

        counter += 1
        if quartet.issubset(frozenset_of_taxa):  # if the tree contains the quartet
            p = round((counter / len(quartet_dictionary)) * 100, 2)
            e = round(extraction_needed/len(quartet_dictionary) * 100, 2)
            if timing:
                sys.stdout.write("                         Tree support progress: %f%% \t Extractions Needed: %f%%   \r" % (p, e) )
                sys.stdout.flush()
            try:
                dict_index = quartetBipartitionSupportHelper(tree, quartet_dictionary, quartet, bipartition_encoding, taxon_label_list, bitstring_encoding, bipartition_dictionary)
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

    dict_index = manualBitmaskSearchV2(sorted_quartet, bipartition_dictionary)
    if (dict_index >= 0):
        quartet_dictionary[quartet][dict_index] = quartet_dictionary[quartet][dict_index] + 1
        return dict_index


    # ERROR
    if dict_index < 0:
        return dict_index
    raise Exception('Error: Topology is not a match')

def manualBitmaskSearchV2(sorted_quartet, bipartition_dictionary):

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[1]])
    second_duet = frozenset([sorted_quartet[2], sorted_quartet[3]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 0
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 0

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[2]])
    second_duet = frozenset([sorted_quartet[1], sorted_quartet[3]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 1
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 1

    first_duet  = frozenset([sorted_quartet[0], sorted_quartet[3]])
    second_duet = frozenset([sorted_quartet[1], sorted_quartet[2]])
    if len(bipartition_dictionary[first_duet]['1'].intersection(bipartition_dictionary[second_duet]['0'])) is not 0:
        return 2
    if len(bipartition_dictionary[second_duet]['1'].intersection(bipartition_dictionary[first_duet]['0'])) is not 0:
        return 2

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

def quartetExtractionSupportHelper(tree, quartet_dictionary, quartet):
    sorted_quartet = list(quartet)
    sorted_quartet.sort()

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

    raise Exception('Error: Topology is not a match')


# Input: A Dendropy TreeList
# Output: Writes the output quartet_dictionary to a file
def writeQuartetDictionaries(bootstrap_tree_list, output_filepath, verbose=False, quiet=False, timing=False):
    if not quiet:
        print("Writing out gene tree data...")
        print()

    start = time.perf_counter()
    if timing:
        print('START: ', start)

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

    file_obj = open(output_filepath, 'wb')
    pickle.dump(quartet_dictionary, file_obj)
    file_obj.close()


def runProgram(sampleTree, output_directory, verbose=False, quiet=False, timing=False):
    output_filename = sampleTree.split('/')[-1]
    output_filepath = output_directory + output_filename + '.quartet_dictionary'
    if verbose:
        print()
        print("Sample Tree: ", sampleTree)
        print("Output File: ", output_filepath)

    if timing:
        verbose = False

    writeQuartetDictionaries(readTree(sampleTree, quiet), output_filepath, verbose, quiet, timing)


    # file_obj = open(output_filepath, 'rb')
    # quartet_dictionary = pickle.load(file_obj)
    # print()
    # print("Full quartet dictionary from pickle:")
    # [print(quartet, quartet_dictionary[quartet]) for quartet in quartet_dictionary]
    # print()
    # file_obj.close()


runProgram('test_trees/highest_support.txt', 'quartet_dictionaries/', verbose=True, timing=True)
runProgram('test_trees/medium_support.txt', 'quartet_dictionaries/', verbose=True, timing=True)
runProgram('test_trees/low_support.txt', 'quartet_dictionaries/', verbose=True, timing=True)
# runProgram('run_files/RAxML_bootstrap.orfg1.last_2.subSample', 'quartet_dictionaries/', verbose=True, timing=True)
