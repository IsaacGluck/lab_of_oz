# Evolutionary Computational Genomics Lab
## Isaac Gluck '19, Junior Research Scholar for Professor Zhaxybayeva

Python Scripts; working towards a *Novel Algorithm to Detect Exchange of Genetic Material Among Microorganisms*.

### Required Libraries
[ETE Toolkit](https://www.etetoolkit.org)
[Dendropy](https://dendropy.org/)

### Pseudocode

1. Read in bootstrap sample gene trees from a file, return dendropy trees
```python
def readTreeFile(filename):
    file = open_file(filename)
    trees = []
    for line in file:
      trees.add(line)
    return trees
```

2. Take a list of trees, return all unique combinations of 4 taxa as keys (using tuples) in a dictionary
  - **Dictionary Structure** { (a, b, c, d): [t1, b1, t2, b2, t3, b3] }
    - a-d are taxa names, t1-3 are the possible topologies and b1-3 are their bootstrap support values.

```python
def getQuartets(trees):
  quartet_dictionary = {}

  for tree in trees:
    list_of_quartets = get_all_combinations_of_four(tree)
    for quartet in list_of_quartets:
      quartet_dictionary.add(quartet, get_topologies(quartet))

  return quartet_dictionary
```
