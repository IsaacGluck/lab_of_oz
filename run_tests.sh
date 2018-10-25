#!/bin/bash

./MethodsV2.py test_trees/reference_tree.txt test_trees/highest_support.txt test_trees/highest_support.txt -v -c 8 -o test_results/highest_support_tree.tre > test_results/highest_support_output.txt

./MethodsV2.py test_trees/reference_tree.txt test_trees/high_conflict.txt test_trees/high_conflict.txt -v -c 8 -o test_results/high_conflict_tree.tre > test_results/high_conflict_output.txt

./MethodsV2.py test_trees/reference_tree.txt test_trees/low_support.txt test_trees/low_support.txt -v -c 8 -o test_results/low_support_tree.tre > test_results/low_support_output.txt

./MethodsV2.py test_trees/reference_tree.txt test_trees/medium_support.txt test_trees/medium_support.txt -v -c 8 -o test_results/medium_support_tree.tre > test_results/medium_support_output.txt
