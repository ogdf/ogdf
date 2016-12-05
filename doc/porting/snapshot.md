[OGDF](README.md) » [Porting Guide](porting.md) » Snapshot

# Porting from Baobab to current snapshot

## General

### C++11

OGDF now requires C++11 features.
Make sure you do not use another standard when compiling your user programs.


### Compiling for debug mode

If you want to compile for debug mode (even if it is user code, linking
against OGDF), you do not need to define `OGDF_DEBUG` any longer. This is
automatically defined using CMake.


### Header files

It might be necessary that you need to include additional header files
for code that worked before. For example, `#include <ogdf/basic/NodeArray.h>`
does not ensure that you also have an `EdgeArray`.

We also removed header files that contained declarations of classes whose
implementation has already been removed in former releases.
This includes `FixedUpwardEmbeddingInserter.h` and `UpwardPlanarizationLayout.h`.

### Global namespace

In former versions of the OGDF some symbols were added to the global namespace.
This includes but isn't limited to `ifstream`, `ofstream`, `min`, `max`, and `numeric_limits`.
Thus, you may be required to explicitly use the OGDF namesapce in some places.


### Exceptions

When a `PreconditionViolationException` or `AlgorithmFailureException` was only thrown
in debug mode, it is now replaced by an assertion (`OGDF_ASSERT`).


### COIN-OR

We have removed Symphony and Cgl from the included COIN-OR package.
The only LP solver we ship is hence Clp.

If your code used Cgl or Symphony, you have to link your objects to
the original versions of these libraries now.


## Macros

Some macros have been removed because they are not necessary any longer.

This includes
 * `OGDF_NEW` (simply use `new` instead)
 * `ABACUS_LP_OSI` (Abacus always uses COIN)
 * `OGDF_NO_COMPILER_TLS` (thread-local storage is a C++11 feature)
 * `OGDF_DECL_THREAD` (simply use `thread_local` instead)

We also have removed most of the iteration macros since C++11 offers range-based for loops.
This includes the removal of `forall_nodes` and `forall_edges`.

```cpp
node v;
Graph G;

forall_nodes(v, G) { ... }
```
should be replaced by
```cpp
for(node v : G.nodes) { ... }
```

There must now be a semicolon after usage of the macro `OGDF_ASSERT` and `OGDF_ASSERT_IF`.

## Changed class, function, header names

### stNumber(), testSTnumber()

To access the st-numbering functions, include `ogdf/basic/STNumbering.h`
(instead of `ogdf/basic/extended_graph_alg.h`).

The respective functions were renamed:

| Former       | New                 |
|--------------|---------------------|
| stNumber     | computeSTNumbering  |
| testSTnumber | isSTNumbering       |

## Graph class

### List of Adjacent Edges

`NodeElement::adjEdges` was renamed to `NodeElement::adjEntries`.
This means you have to change `for(adjEntry adj : v->adjEdges)` to `for(adjEntry adj : v->adjEntries)`.

The following getter-methods were moved from `Graph` to `NodeElement`.
All of these methods take a single list as the only parameter now.
 * `adjEntries()` (renamed to `allAdjEntries()`)
 * `adjEdges()`
 * `inEdges()`
 * `outEdges()`


### Hiding and restoring edges

Hiding and restoring edges was not a safe operation in former releases.
Replace code like
```c++
    G.hideEdge(e);
    // ...
    G.restoreEdge(e);
```
by
```c++
   Graph::HiddenEdgeSet hiddenSet(G);
   hiddenSet.hide(e);
   // ...
   hiddenSet.restore(e);
```
All edges are restored by `hiddenSet.restore()` or automatically by the `HiddenEdgeSet` destructor.


### Typofix: collaps

The method `collaps` has been renamed to `collapse`.


## GraphCopy class

The `newEdge()` methods to add edges at predefined positions in the adjacency list
have been removed from `GraphCopy`.
You can instead use the respective `newEdge()` function inherited from `Graph`,
followed by `GraphCopy::setEdge()`.


## List & ListPure class

The deprecated methods `List::exchange()` and `ListPure::exchange()` were removed.
Use `List::swap()` and `ListPure::swap()` instead.

## Planar Subgraph Algorithms

`BoyerMyrvoldSubgraph` was renamed to `PlanarSubgraphBoyerMyrvold` and `FastPlanarSubgraph` to `PlanarSubgraphFast`.
`MaximalPlanarSubgraphSimple` now supports arbitrary start heuristics.
As such, the `doMaximize` option was removed from `PlanarSubgraphBoyerMyrvold` (formerly `BoyerMyrvoldSubgraph`).
`MaximalPlanarSubgraphSimple::clone()` should be used to obtain a copy of the respective instance since `MaximalPlanarSubgraphSimple(const MaximalPlanarSubgraphSimple &mps)` is no longer available.

Support for edge weights was added to `MaximumCPlanarSubgraph`. The signature of `call` includes the new parameter `EdgeArray<int> *pCost` that can be set to `nullptr` for uniform weight.

The Module `PlanarSubgraphModule` is now a template. The implementations `PlanarSubgraphCactus`, `MaximalPlanarSubgraphSimple` and `PlanarSubgraphEmpty` are templates itself. All other implementations of the module inherit from `PlanarSubgraphModule<int>`.

### MaximumCPlanarSubgraph

`MaximumCPlanarSubgraph` supports advanced calls that also return the edges required to connect the input graph.
To avoid conflicts with the method defined by `CPlanarSubgraphModule`, `MaximumCPlanarSubgraph::call` was renamed to
`MaximumCPlanarSubgraph::callAndConnect`.

## Heaps

In an attempt to unify the existing heap classes, `BinaryHeap` and `BinaryHeap2` were replaced by a new `BinaryHeap`.
All heap implementations are derived from the common interface `HeapBase`.
Additionally, the classes `MinPriorityQueue` and `PQueue` were removed.

The following heap implementations were introduced:
 * `BinaryHeap`
 * `BinomialHeap`
 * `FibonacciHeap`
 * `PairingHeap`
 * `RadixHeap`

Priority queues might be realized using the newly introduced `PriorityQueue` idependently of the desired heap implementation.
In contrast to `PriorityQueue` that merely stores directly comparable elements, `PrioritizedQueue` is used to store elements with priorities assigned upon insertion.
For even simpler usage see `PrioritizedMapQueue` that keeps track of handles for the inserted elements but requires elements to be unique.

The tests in `test/src/basic/heap_test.cpp` show exemplary usage of the new classes.

### Method equivalents

#### BinaryHeap
`BinaryHeap` was replaced by `PrioritizedQueue`.
Accessing elements at arbitrary positions is no longer supported.

| Former                    | New                                                      |
|---------------------------|----------------------------------------------------------|
| `BinaryHeap::clear`       | `PrioritizedQueue::clear`                                |
| `BinaryHeap::decPriority` | `PrioritizedQueue::decrease`                             |
| `BinaryHeap::empty`       | `PrioritizedQueue::empty`                                |
| `BinaryHeap::extractMin`  | `PrioritizedQueue::topElement` & `PrioritizedQueue::pop` |
| `BinaryHeap::getMin`      | `PrioritizedQueue::topElement`                           |
| `BinaryHeap::insert`      | `PrioritizedQueue::push`                                 |
| `BinaryHeap::operator[]`  | -                                                        |
| `BinaryHeap::pop`         | `PrioritizedQueue::pop`                                  |
| `BinaryHeap::push`        | `PrioritizedQueue::push`                                 |
| `BinaryHeap::size`        | `PrioritizedQueue::size`                                 |
| `BinaryHeap::top`         | `PrioritizedQueue::topElement`                           |

#### BinaryHeap2

`BinaryHeap2` was replaced by `PrioritizedMapQueue`. Querying the capacity of the heap is no longer supported.

| Former                       | New                                                            |
|------------------------------|----------------------------------------------------------------|
| `BinaryHeap2::capacity`      | -                                                              |
| `BinaryHeap2::clear`         | `PrioritizedMapQueue::clear`                                   |
| `BinaryHeap2::decreaseKey`   | `PrioritizedMapQueue::decrease`                                |
| `BinaryHeap2::empty`         | `PrioritizedMapQueue::empty`                                   |
| `BinaryHeap2::extractMin`    | `PrioritizedMapQueue::topElement` & `PrioritizedMapQueue::pop` |
| `BinaryHeap2::getMinimumKey` | `PrioritizedMapQueue::topPriority`                             |
| `BinaryHeap2::getPriority`   | `PrioritizedMapQueue::priority`                                |
| `BinaryHeap2::insert`        | `PrioritizedMapQueue::push`                                    |
| `BinaryHeap2::size`          | `PrioritizedMapQueue::size`                                    |
| `BinaryHeap2::topElement`    | `PrioritizedMapQueue::topElement`                              |

#### PQueue

`PQueue` was replaced by `PrioritizedQueue`.

| Former                       | New            |
|------------------------------|----------------|
| `PQueue::del_min`  | `PrioritizedQueue::pop`  |
| `PQueue::find_min` | `PrioritizedQueue::top`  |
| `PQueue::insert`   | `PrioritizedQueue::push` |

#### MinPriorityQueue

`MinPriorityQueue` was replaced by `PrioritizedQueue`. Note that `PrioritizedQueue::size` returns the number of elements in the heap instead of the capacity.

| Former                               | New                          |
|--------------------------------------|------------------------------|
| `MinPriorityQueue::count`            | `PrioritizedQueue::size`     |
| `MinPriorityQueue::decreasePriority` | `PrioritizedQueue::decrease` |
| `MinPriorityQueue::empty`            | `PrioritizedQueue::empty`    |
| `MinPriorityQueue::getMin`           | `PrioritizedQueue::top`      |
| `MinPriorityQueue::insert`           | `PrioritizedQueue::push`     |
| `MinPriorityQueue::pop`              | `PrioritizedQueue::pop`      |
| `MinPriorityQueue::size`             | -                            |

## Layout algorithms and graph constraints

Layout algorithms no longer support `GraphConstraints`.
Use `call(GraphAttributes)` instead of `call(GraphAttributes, GraphConstraints)`.
Graph constraints were removed entirely.
Removed classes: `Constraint`, `ConstraintManager`, and `GraphConstraints`.

## GraphIO

Methods that take filenames (instead of streams) are deprecated and will be removed in future versions.
For example, instead of using

```c++
	GraphIO::readLEDA(G, "circulant.lgr");
	GraphIO::writeChaco(G, "circulant.gml");
```

you have to create the streams manually, i.e.,

```c++
	ifstream is("circulant.lgr");
	GraphIO::readLEDA(G, is);
	ofstream os("circulant.gml");
	GraphIO::writeGML(G, os);
```

There are also helper functions `GraphIO::read()` and `GraphIO::write()`
that accept filenames for `GraphAttribute`s and `Graph`s,
so you can also write

```c++
	GraphIO::read(G, "circulant.lgr", GraphIO::readLEDA);
	GraphIO::write(G, "circulant.gml", GraphIO::writeGML);
```

Do not confuse the new filename-based `GraphIO::read()` helper function with
our new (experimental) generic reader function of the same name.
However, this generic reader even allows to write

```c++
	GraphIO::read(G, "circulant.lgr");
	GraphIO::write(G, "circulant.gml", GraphIO::writeGML);
```

## Filesystem functions

If you have used filesystem functions in your code,
you can include `ogdf/basic/filesystem.h` now and you get them back.
They are marked deprecated now since there are maintained and portable
alternatives like [tinydir](https://github.com/cxong/tinydir).
Also C++17 will have filesystem functions.

## CombinatorialEmbedding

`CombinatorialEmbedding::splitFace(adjEntry, node)` and `CombinatorialEmbedding::splitFace(node, adjEntry)`
are used to insert degree-0 nodes without having to call `CombinatorialEmbedding::computeFaces()` again.
In previous version both methods could be used with non-isolated nodes
if the last adjacency entry belonged to the same face as the other adjacency entry.
This is no longer supported and the node has to be isolated.
To reflect this functional change the respective versions of `splitFace` were renamed to `addEdgeToIsolatedNode`.
Use `myEmbedding.splitFace(v->lastAdj(), adj)` instead of `myEmbedding.splitFace(v, adj)` if `v` isn't isolated (and you are sure that the adjacency entries lie on the same face!).
Use `myEmbedding.addEdgeToIsolatedNode(v, adj)` if `v` is isolated.

## OGDF_USE_SSE2, OGDF_CHECK_SSE2

The macros mentioned above have been removed.
You can check for the SSE2 CPU feature directly using
`ogdf::System::cpuSupports(ogdf::cpufSSE2)`.

## ModuleOption

The template class `ModuleOption` was removed. `std::unique_ptr` should be used instead.

## Angle functions for DPoint and GenericPoint

The static methods `angle` and `angleDegrees` of `DPoint` and `GenericPoint` are now
non-static. That means, old calls of

```c++
	double angle = DPoint::angle(point1, point2, point3);
```

now become

```c++
	double angle = point1.angle(point2, point3);
```
