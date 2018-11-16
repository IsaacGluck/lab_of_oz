#!/usr/bin/env python
from dendropy import Tree, TreeList
import re

def addSiblingsFromLabels(concat_tree_file, species_tree_file, split_string):
    concat_tree = Tree.get(path=concat_tree_file, schema="newick", preserve_underscores=True)
    species_tree = Tree.get(path=species_tree_file, schema="newick", preserve_underscores=True)

    concat_tree_leaves = [leaf for leaf in concat_tree.leaf_nodes()]

    for leaf_node in species_tree.leaf_nodes():
        species_label = leaf_node.taxon.label
        # print([c.taxon for c in leaf_node.sibling_nodes()])
        similar_nodes = find_similar(concat_tree_leaves, species_label)
        if (len(similar_nodes) is 1):
            string_to_add = split_string + similar_nodes[0].taxon.label.split(split_string)[1]
            new_label = leaf_node.taxon.label + string_to_add
            leaf_node.taxon.label = new_label
            continue
            # print('old: ', leaf_node.taxon.label, '   new: ', new_label)
        elif (len(similar_nodes) > 1):
            # print(leaf_node._parent_node._child_nodes)
            parent = leaf_node._parent_node

            for new_node in similar_nodes:
                parent.add_child(new_node)
            old_node = parent.remove_child(leaf_node)
            # print('old node: ', old_node, '       new nodes: ', parent.child_nodes())
            # print()
            continue

    print(species_tree.as_string(schema="newick", suppress_internal_node_labels=False))




def find_similar(leaves, species_label):
    similar = []

    p = re.compile(species_label)
    for leaf in leaves:
        if p.match(leaf.taxon.label) is not None:
            similar.append(leaf)

    return similar


addSiblingsFromLabels('./run_files/RAxML_bestTree.rcGTA_cat', './reference_tree/reference_tree_LCsubset.tre', '_GTAclstr_')
