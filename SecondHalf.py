#!/usr/bin/env python
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from itertools import combinations
from math import log
from pprint import pprint
import pickle, sys, os, time, argparse


# Input: list of files with serialized quartet_dictionaries
# Output: A list of quartet_dictionary objects
def readPickledTrees(quartetDictionaryFileList, quiet, timing):
    if not quiet:
        print()
        print("Reading in quartet dictionary files...")
        print()

    quartet_dictionary_list = []

    counter = 0
    for filename in quartetDictionaryFileList:

        file_obj = open(filename, 'rb')
        qd = pickle.load(file_obj)
        quartet_dictionary_list.append(qd)
        file_obj.close()

        counter += 1
        p = round((counter / len(quartetDictionaryFileList)) * 100, 2)
        if timing:
            sys.stdout.write("[%d/%d] Unpickling file progress: %f%%   \r" % (counter, len(quartetDictionaryFileList), p) )
            sys.stdout.flush()

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

def buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree="output_tree.tre", quiet=False, timing=False):
    reference_tree = Tree.get(path=referenceTreeFile, schema="newick", preserve_underscores=True)
    reference_tree.is_rooted = True
    reference_tree.encode_bipartitions()
    bipartition_encoding = set(b.split_bitmask for b in reference_tree.bipartition_encoding)
    bitstring_encoding = []
    for b in reference_tree.bipartition_encoding:
        if not b.is_trivial():
            bitstring_encoding.append(b.split_as_bitstring())
    taxon_label_list = [(n.taxon.label) for n in reference_tree.leaf_nodes()]
    bipartition_dictionary = makeBipartitionDictionary(taxon_label_list, bitstring_encoding)




    tn = reference_tree.taxon_namespace

    start = start = time.perf_counter()
    if timing:
        print('START BUILD LABEL TREE: ', start)
    splits = getListOfSplits(tn, reference_tree)
    if timing:
        print('GOT LIST OF SPLITS: ', (time.perf_counter() - start))


    total_exist = 1
    not_exist = 0
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
                    sys.stdout.write("Build Label Tree Progress [%s/%s] : %s%%   %s%%          not found: %f%  %d \r" % (str(counter), str(len(splits)), p, inner_p, not_exist/total_exist * 100, not_exist))
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

                    # Backup
                    if ((not results[0]) and (not results[1]) and (not results[2])):
                        result_index = manualBitmaskSearchV2(combined_taxa_labels, bipartition_dictionary)
                        results[result_index] = True

                    total_exist += 1
                    if ((not results[0]) and (not results[1]) and (not results[2])):
                        not_exist += 1

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
        print("Quartet Dictionary File List: ", quartetDictionaryFileList)
        print("Bootstrap Cutoff Value: ", bootstrap_cutoff_value)
        print("Output Tree File: ", output_tree)

    if timing:
        verbose = False

    try:
        reference_tree = Tree.get(path=referenceTreeFile, schema="newick", preserve_underscores=True)
    except:
        print("Error with file '{}': please only use files with newick tree format".format(referenceTreeFile))
        sys.exit()

    quartet_dictionary_list = readPickledTrees(quartetDictionaryFileList, quiet, timing)

    full_quartet_dictionary = buildFullSupport(quartet_dictionary_list, bootstrap_cutoff_value, verbose, quiet, timing)
    if verbose:
        print("Full quartet dictionary with support values")
        [print(quartet, full_quartet_dictionary[quartet])
               for quartet in full_quartet_dictionary]
        print()
    buildLabeledTree(referenceTreeFile, full_quartet_dictionary, output_tree, quiet, timing)

# quartetDictionaryFileList = ['quartet_dictionaries/highest_support.txt.quartet_dictionary', 'quartet_dictionaries/low_support.txt.quartet_dictionary', 'quartet_dictionaries/medium_support.txt.quartet_dictionary']
# ./MethodsV2.py test_trees/reference_tree.txt test_trees/highest_support.txt test_trees/low_support.txt test_trees/medium_support.txt -t -c 8
# runProgram('test_trees/reference_tree.txt', quartetDictionaryFileList, bootstrap_cutoff_value=8, timing=True)



parser = argparse.ArgumentParser()
parser.add_argument("reference_tree_file", metavar='<Reference Tree File>', help="The path of the reference tree file")
parser.add_argument('quartet_dictionary_file_list', metavar='<Quartet Dictionary Files>', nargs='+',
                    help='The gene tree file paths containing bootstrap trees')
parser.add_argument("-v", "--verbose", action="store_true", default=False)
parser.add_argument("-q", "--quiet", action="store_true", default=False)
parser.add_argument("-t", "--timing", action="store_true", default=False, help='Setting timing turns verbose off.')
parser.add_argument("-c", "--cutoff", default=80, type=int,
                    help="Bootstrap Cutoff Value")
parser.add_argument("-o", "--output_file", default="output_tree.tre",
                    help="Output file for resulting tree with support. The default is 'output_tree.tre'")
args = parser.parse_args()

runProgram(args.reference_tree_file, args.quartet_dictionary_file_list,
           bootstrap_cutoff_value=args.cutoff, output_tree=args.output_file, verbose=args.verbose, quiet=args.quiet, timing=args.timing)



# REGEX TO REMOVE BRANCH LENGTHS
# :\d+\.\d+(e-\d+)?

# QUICK RUN
# time ./SecondHalf.py run_files/RAxML_bestTree.rcGTA_cat quartet_dictionaries/RAxML_bootstrap.orfg1.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg10_5.last_3.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg3_5.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg7.last_2.subSample.quartet_dictionary -c 8 -t -v

# FULL RUN
# time ./SecondHalf.py run_files/RAxML_bestTree.rcGTA_cat quartet_dictionaries/RAxML_bootstrap.orfg1.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg10.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg10_5.last_3.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg11.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg12.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg14.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg15.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg2.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg3.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg3_5.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg4.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg5.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg6.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg7.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg8.last_2.subSample.quartet_dictionary quartet_dictionaries/RAxML_bootstrap.orfg9.last_2.subSample.quartet_dictionary  -c 8 -t -v
