# Reference Tree:
```
(((A,B),E),(D,C));
```

### Highest Support:

```
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(D,C));
```

### Run Output Highest Support, cutoff=8
```
Reading in files...
Combining gene tree data into one dictionary...
Making the quartet dictionary...

Full quartet dictionary:
frozenset({'E', 'B', 'A', 'D'}) [<Tree object at 0x103198be0>, 10, <Tree object at 0x103198d68>, 0, <Tree object at 0x1031989b0>, 0]
frozenset({'C', 'E', 'B', 'A'}) [<Tree object at 0x103198e80>, 10, <Tree object at 0x1031989e8>, 0, <Tree object at 0x103198b00>, 0]
frozenset({'C', 'B', 'A', 'D'}) [<Tree object at 0x10318c1d0>, 10, <Tree object at 0x1031a8b38>, 0, <Tree object at 0x103198da0>, 0]
frozenset({'C', 'E', 'A', 'D'}) [<Tree object at 0x10318c208>, 0, <Tree object at 0x1031a8e80>, 0, <Tree object at 0x1031ae898>, 10]
frozenset({'C', 'E', 'B', 'D'}) [<Tree object at 0x10318c048>, 0, <Tree object at 0x103198c18>, 0, <Tree object at 0x1031ae208>, 10]

Making the quartet dictionary...

Full quartet dictionary:
frozenset({'E', 'B', 'A', 'D'}) [<Tree object at 0x1031bb710>, 10, <Tree object at 0x1031bbc18>, 0, <Tree object at 0x1031bb9e8>, 0]
frozenset({'C', 'E', 'B', 'A'}) [<Tree object at 0x1031bb7b8>, 10, <Tree object at 0x1031c5048>, 0, <Tree object at 0x1031bb748>, 0]
frozenset({'C', 'B', 'A', 'D'}) [<Tree object at 0x1031bbef0>, 10, <Tree object at 0x1031bbb70>, 0, <Tree object at 0x1031bbcf8>, 0]
frozenset({'C', 'E', 'A', 'D'}) [<Tree object at 0x1031bbb00>, 0, <Tree object at 0x1031bb860>, 0, <Tree object at 0x1031c5e10>, 10]
frozenset({'C', 'E', 'B', 'D'}) [<Tree object at 0x1031c5588>, 0, <Tree object at 0x1031d8828>, 0, <Tree object at 0x1031bb7f0>, 10]

Full quartet dictionary with support values
frozenset({'E', 'B', 'A', 'D'}) [<Tree object at 0x103198be0>, 1.0, <Tree object at 0x103198d68>, 0.0, <Tree object at 0x1031989b0>, 0.0, 1.0]
frozenset({'C', 'E', 'B', 'A'}) [<Tree object at 0x103198e80>, 1.0, <Tree object at 0x1031989e8>, 0.0, <Tree object at 0x103198b00>, 0.0, 1.0]
frozenset({'C', 'B', 'A', 'D'}) [<Tree object at 0x10318c1d0>, 1.0, <Tree object at 0x1031a8b38>, 0.0, <Tree object at 0x103198da0>, 0.0, 1.0]
frozenset({'C', 'E', 'A', 'D'}) [<Tree object at 0x10318c208>, 0.0, <Tree object at 0x1031a8e80>, 0.0, <Tree object at 0x1031ae898>, 1.0, 1.0]
frozenset({'C', 'E', 'B', 'D'}) [<Tree object at 0x10318c048>, 0.0, <Tree object at 0x103198c18>, 0.0, <Tree object at 0x1031ae208>, 1.0, 1.0]

{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0ac8>}
{'right': [<Taxon 0x1031d0da0 'D'>, <Taxon 0x1031bbef0 'C'>], 'left': [<Taxon 0x1031c5588 'A'>, <Taxon 0x1031c57b8 'B'>, <Taxon 0x1031d05f8 'E'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0ef0>}
{'right': [<Taxon 0x1031d05f8 'E'>, <Taxon 0x1031d0da0 'D'>, <Taxon 0x1031bbef0 'C'>], 'left': [<Taxon 0x1031c5588 'A'>, <Taxon 0x1031c57b8 'B'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0f60>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0cf8>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0cc0>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d0630>}
{'right': [<Taxon 0x1031c5588 'A'>, <Taxon 0x1031c57b8 'B'>, <Taxon 0x1031d05f8 'E'>], 'left': [<Taxon 0x1031d0da0 'D'>, <Taxon 0x1031bbef0 'C'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031c5128>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031c5e48>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1031d45c0>}

[&R] (((A,B)1.0,E)1.0,(D,C)1.0);
```


# Lowest Support:
```
(((A,B),E),(D,C));
(((A,C),E),(D,B));
(((A,D),E),(B,C));
(((E,B),D),(A,C));
(((A,C),B),(D,E));
(((E,C),A),(D,B));
(((E,D),C),(A,B));
(((C,A),B),(E,D));
(((D,C),A),(B,E));
(((E,D),C),(B,A));
```

### Run Output Lowest Support, cutoff=8
```
Reading in files...
Combining gene tree data into one dictionary...
Making the quartet dictionary...

Full quartet dictionary:
frozenset({'B', 'D', 'E', 'A'}) [<Tree object at 0x1035a6be0>, 5, <Tree object at 0x1035a6d68>, 3, <Tree object at 0x1035a69b0>, 2]
frozenset({'C', 'B', 'E', 'A'}) [<Tree object at 0x1035a6e80>, 4, <Tree object at 0x1035a69e8>, 5, <Tree object at 0x1035a6b00>, 1]
frozenset({'C', 'B', 'D', 'A'}) [<Tree object at 0x10359a1d0>, 4, <Tree object at 0x1035b6b38>, 5, <Tree object at 0x1035a6da0>, 1]
frozenset({'C', 'D', 'E', 'A'}) [<Tree object at 0x10359a208>, 6, <Tree object at 0x1035b6e80>, 2, <Tree object at 0x1035bc898>, 2]
frozenset({'C', 'B', 'D', 'E'}) [<Tree object at 0x10359a048>, 5, <Tree object at 0x1035a6c18>, 2, <Tree object at 0x1035bc208>, 3]

Making the quartet dictionary...

Full quartet dictionary:
frozenset({'B', 'D', 'E', 'A'}) [<Tree object at 0x1035c9860>, 5, <Tree object at 0x1035c9da0>, 3, <Tree object at 0x1035c96a0>, 2]
frozenset({'C', 'B', 'E', 'A'}) [<Tree object at 0x1035c9828>, 4, <Tree object at 0x1035d31d0>, 5, <Tree object at 0x1035c9e10>, 1]
frozenset({'C', 'B', 'D', 'A'}) [<Tree object at 0x1035c9be0>, 4, <Tree object at 0x1035c9f28>, 5, <Tree object at 0x1035c9b38>, 1]
frozenset({'C', 'D', 'E', 'A'}) [<Tree object at 0x1035c9b00>, 6, <Tree object at 0x1035c9b70>, 2, <Tree object at 0x1035d3208>, 2]
frozenset({'C', 'B', 'D', 'E'}) [<Tree object at 0x1035d35c0>, 5, <Tree object at 0x1035e6588>, 2, <Tree object at 0x1035c99e8>, 3]

Full quartet dictionary with support values
frozenset({'B', 'D', 'E', 'A'}) [<Tree object at 0x1035a6be0>, 0.0, <Tree object at 0x1035a6d68>, 0.0, <Tree object at 0x1035a69b0>, 0.0, 1.0]
frozenset({'C', 'B', 'E', 'A'}) [<Tree object at 0x1035a6e80>, 0.0, <Tree object at 0x1035a69e8>, 0.0, <Tree object at 0x1035a6b00>, 0.0, 1.0]
frozenset({'C', 'B', 'D', 'A'}) [<Tree object at 0x10359a1d0>, 0.0, <Tree object at 0x1035b6b38>, 0.0, <Tree object at 0x1035a6da0>, 0.0, 1.0]
frozenset({'C', 'D', 'E', 'A'}) [<Tree object at 0x10359a208>, 0.0, <Tree object at 0x1035b6e80>, 0.0, <Tree object at 0x1035bc898>, 0.0, 1.0]
frozenset({'C', 'B', 'D', 'E'}) [<Tree object at 0x10359a048>, 0.0, <Tree object at 0x1035a6c18>, 0.0, <Tree object at 0x1035bc208>, 0.0, 1.0]

{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035d3208>}
{'right': [<Taxon 0x1035ded30 'D'>, <Taxon 0x1035c9f28 'C'>], 'left': [<Taxon 0x1035de7f0 'A'>, <Taxon 0x1035dec88 'B'>, <Taxon 0x1035de4e0 'E'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035d37b8>}
{'right': [<Taxon 0x1035de4e0 'E'>, <Taxon 0x1035ded30 'D'>, <Taxon 0x1035c9f28 'C'>], 'left': [<Taxon 0x1035de7f0 'A'>, <Taxon 0x1035dec88 'B'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035dea58>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035de748>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035de940>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035de5c0>}
{'right': [<Taxon 0x1035de7f0 'A'>, <Taxon 0x1035dec88 'B'>, <Taxon 0x1035de4e0 'E'>], 'left': [<Taxon 0x1035ded30 'D'>, <Taxon 0x1035c9f28 'C'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035deb00>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035de438>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1035e2518>}

[&R] (((A,B)'-1.0',E)'-1.0',(D,C)'-1.0');
```

# High Conflict:
```
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
(((D,A),C),(E,B));
```

### Run Output High Conflict, cutoff=8
```
Reading in files...
Combining gene tree data into one dictionary...
Making the quartet dictionary...

Full quartet dictionary:
frozenset({'D', 'E', 'C', 'A'}) [<Tree object at 0x10e3f4be0>, 0, <Tree object at 0x10e3f4d68>, 10, <Tree object at 0x10e3f49b0>, 0]
frozenset({'D', 'B', 'C', 'A'}) [<Tree object at 0x10e3f4e80>, 0, <Tree object at 0x10e3f49e8>, 0, <Tree object at 0x10e3f4b00>, 10]
frozenset({'D', 'B', 'E', 'A'}) [<Tree object at 0x10e3e81d0>, 0, <Tree object at 0x10e404b38>, 10, <Tree object at 0x10e3f4da0>, 0]
frozenset({'D', 'E', 'B', 'C'}) [<Tree object at 0x10e3e8208>, 0, <Tree object at 0x10e404e80>, 0, <Tree object at 0x10e40a898>, 10]
frozenset({'E', 'B', 'C', 'A'}) [<Tree object at 0x10e3e8048>, 0, <Tree object at 0x10e3f4c18>, 10, <Tree object at 0x10e40a208>, 0]

Making the quartet dictionary...

Full quartet dictionary:
frozenset({'D', 'E', 'C', 'A'}) [<Tree object at 0x10e417a90>, 0, <Tree object at 0x10e417b38>, 10, <Tree object at 0x10e417ac8>, 0]
frozenset({'D', 'B', 'C', 'A'}) [<Tree object at 0x10e417ef0>, 0, <Tree object at 0x10e4219e8>, 0, <Tree object at 0x10e4177b8>, 10]
frozenset({'D', 'B', 'E', 'A'}) [<Tree object at 0x10e417b70>, 0, <Tree object at 0x10e417ba8>, 10, <Tree object at 0x10e417b00>, 0]
frozenset({'D', 'E', 'B', 'C'}) [<Tree object at 0x10e417f98>, 0, <Tree object at 0x10e4176a0>, 0, <Tree object at 0x10e421cf8>, 10]
frozenset({'E', 'B', 'C', 'A'}) [<Tree object at 0x10e421b38>, 0, <Tree object at 0x10e434908>, 10, <Tree object at 0x10e417908>, 0]

Full quartet dictionary with support values
frozenset({'D', 'E', 'C', 'A'}) [<Tree object at 0x10e3f4be0>, 0.0, <Tree object at 0x10e3f4d68>, 1.0, <Tree object at 0x10e3f49b0>, 0.0, 1.0]
frozenset({'D', 'B', 'C', 'A'}) [<Tree object at 0x10e3f4e80>, 0.0, <Tree object at 0x10e3f49e8>, 0.0, <Tree object at 0x10e3f4b00>, 1.0, 1.0]
frozenset({'D', 'B', 'E', 'A'}) [<Tree object at 0x10e3e81d0>, 0.0, <Tree object at 0x10e404b38>, 1.0, <Tree object at 0x10e3f4da0>, 0.0, 1.0]
frozenset({'D', 'E', 'B', 'C'}) [<Tree object at 0x10e3e8208>, 0.0, <Tree object at 0x10e404e80>, 0.0, <Tree object at 0x10e40a898>, 1.0, 1.0]
frozenset({'E', 'B', 'C', 'A'}) [<Tree object at 0x10e3e8048>, 0.0, <Tree object at 0x10e3f4c18>, 1.0, <Tree object at 0x10e40a208>, 0.0, 1.0]

{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e4219e8>}
{'right': [<Taxon 0x10e42cba8 'D'>, <Taxon 0x10e417b70 'C'>], 'left': [<Taxon 0x10e42c748 'A'>, <Taxon 0x10e42c3c8 'B'>, <Taxon 0x10e42cc50 'E'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e421780>}
{'right': [<Taxon 0x10e42cc50 'E'>, <Taxon 0x10e42cba8 'D'>, <Taxon 0x10e417b70 'C'>], 'left': [<Taxon 0x10e42c748 'A'>, <Taxon 0x10e42c3c8 'B'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e42c940>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e42cb38>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e42ccf8>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e42c4e0>}
{'right': [<Taxon 0x10e42c748 'A'>, <Taxon 0x10e42c3c8 'B'>, <Taxon 0x10e42cc50 'E'>], 'left': [<Taxon 0x10e42cba8 'D'>, <Taxon 0x10e417b70 'C'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e42cbe0>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e3d6a20>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x10e430828>}

[&R] (((A,B)'-1.0',E)'-0.3333333333333333',(D,C)'-0.3333333333333333');
```



# Medium Support/Structure Changes:
```
(((A,B),E),(D,C));
((A,B),((D,C),E));
(((A,B),D),(E,C));
(((A,B),(D,E)),C);
(((A,B),C),(D,E));
(((E,B),A),(D,C));
(((A,B),E),(D,C));
(((A,B),E),(C,D));
((A,E),((D,C),B));
(((A,B),E),(D,C));
```

### Run Output Medium Support, cutoff=8
```
Reading in files...
Combining gene tree data into one dictionary...
Making the quartet dictionary...

Full quartet dictionary:
frozenset({'B', 'E', 'D', 'A'}) [<Tree object at 0x11088fbe0>, 8, <Tree object at 0x11088fd68>, 1, <Tree object at 0x11088f9b0>, 1]
frozenset({'B', 'E', 'A', 'C'}) [<Tree object at 0x11088fe80>, 8, <Tree object at 0x11088f9e8>, 1, <Tree object at 0x11088fb00>, 1]
frozenset({'B', 'D', 'A', 'C'}) [<Tree object at 0x1108831d0>, 10, <Tree object at 0x11089fb38>, 0, <Tree object at 0x11088fda0>, 0]
frozenset({'E', 'D', 'A', 'C'}) [<Tree object at 0x110883208>, 2, <Tree object at 0x11089fe80>, 1, <Tree object at 0x1108a5898>, 7]
frozenset({'B', 'E', 'D', 'C'}) [<Tree object at 0x110883048>, 2, <Tree object at 0x11088fc18>, 1, <Tree object at 0x1108a5208>, 7]

Making the quartet dictionary...

Full quartet dictionary:
frozenset({'B', 'E', 'D', 'A'}) [<Tree object at 0x1108b2358>, 8, <Tree object at 0x1108b2978>, 1, <Tree object at 0x1108b2748>, 1]
frozenset({'B', 'E', 'A', 'C'}) [<Tree object at 0x1108b2ba8>, 8, <Tree object at 0x1108bc3c8>, 1, <Tree object at 0x1108b2cf8>, 1]
frozenset({'B', 'D', 'A', 'C'}) [<Tree object at 0x1108b26a0>, 10, <Tree object at 0x1108b29e8>, 0, <Tree object at 0x1108b2be0>, 0]
frozenset({'E', 'D', 'A', 'C'}) [<Tree object at 0x1108b2b00>, 2, <Tree object at 0x1108b2d30>, 1, <Tree object at 0x1108bce10>, 7]
frozenset({'B', 'E', 'D', 'C'}) [<Tree object at 0x1108bc240>, 2, <Tree object at 0x1108cf828>, 1, <Tree object at 0x1108b2ef0>, 7]

Full quartet dictionary with support values
frozenset({'B', 'E', 'D', 'A'}) [<Tree object at 0x11088fbe0>, 1.0, <Tree object at 0x11088fd68>, 0.0, <Tree object at 0x11088f9b0>, 0.0, 1.0]
frozenset({'B', 'E', 'A', 'C'}) [<Tree object at 0x11088fe80>, 1.0, <Tree object at 0x11088f9e8>, 0.0, <Tree object at 0x11088fb00>, 0.0, 1.0]
frozenset({'B', 'D', 'A', 'C'}) [<Tree object at 0x1108831d0>, 1.0, <Tree object at 0x11089fb38>, 0.0, <Tree object at 0x11088fda0>, 0.0, 1.0]
frozenset({'E', 'D', 'A', 'C'}) [<Tree object at 0x110883208>, 0.0, <Tree object at 0x11089fe80>, 0.0, <Tree object at 0x1108a5898>, 0.0, 1.0]
frozenset({'B', 'E', 'D', 'C'}) [<Tree object at 0x110883048>, 0.0, <Tree object at 0x11088fc18>, 0.0, <Tree object at 0x1108a5208>, 0.0, 1.0]

{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108c6278>}
{'right': [<Taxon 0x1108c6e10 'D'>, <Taxon 0x1108c65c0 'C'>], 'left': [<Taxon 0x1108bc240 'A'>, <Taxon 0x10f6c4470 'B'>, <Taxon 0x1108bcbe0 'E'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108c6ef0>}
{'right': [<Taxon 0x1108bcbe0 'E'>, <Taxon 0x1108c6e10 'D'>, <Taxon 0x1108c65c0 'C'>], 'left': [<Taxon 0x1108bc240 'A'>, <Taxon 0x10f6c4470 'B'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108c6cf8>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108c6d68>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108cb6d8>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108cb860>}
{'right': [<Taxon 0x1108bc240 'A'>, <Taxon 0x10f6c4470 'B'>, <Taxon 0x1108bcbe0 'E'>], 'left': [<Taxon 0x1108c6e10 'D'>, <Taxon 0x1108c65c0 'C'>], 'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108cb5c0>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108cb940>}
{'edge': <dendropy.datamodel.treemodel.Edge object at 0x1108cb668>}

[&R] (((A,B)1.0,E)'-0.3333333333333333',(D,C)'-0.3333333333333333');
```
