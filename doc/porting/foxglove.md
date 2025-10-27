[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Foxglove

# Porting from Elderberry to Foxglove

## Requirements
OGDF now requires C++17 features.
We no longer officially support compilers older than gcc 9, clang 9 or Visual Studio 2017 15.8 (MSVC 19.15).

## Build Process
The file `ogdf/basic/internal/config_autogen.h` has moved from `${PROJECT_BINARY_DIR}/include/`
to `${PROJECT_BINARY_DIR}/include/ogdf-{release,debug}/` (depending on the build type).
If you manually specify include paths for the OGDF make sure to update them accordingly such that
the file can be included as before using `#include <ogdf/basic/internal/config_autogen.h>`.
CMake users should not notice any changes.

Binaries built in debug mode now have the suffix "-debug".

## include-what-you-use (iwyu)
The OGDF now [uses](https://github.com/ogdf/ogdf/pull/239) [iwyu](https://include-what-you-use.org/)
to make sure each source file explicitly lists all header files it uses, but no further unused headers.
This especially means that certain headers that were previously transitively provided but not used by some other header
now might need to be explicitly included in your code. We recommend also running iwyu on your code before (and also after) porting
to get an explicit overview over which imports are used from where.

Due to the `RegisteredArray` changes mentioned below, this also means that the following header files have been removed
as all their functionality is now included in the corresponding `(Cluster/Hyper)Graph.h` file.
```
ogdf/basic/NodeArray.h
ogdf/basic/AdjEntryArray.h
ogdf/basic/FaceArray.h
ogdf/basic/GraphObserver.h
ogdf/cluster/ClusterArray.h
ogdf/cluster/ClusterGraphObserver.h
ogdf/hypergraph/HypergraphArray.h
```
The header `ogdf/basic/NodeSet.h` was replaced by `ogdf/basic/GraphSets.h`, now also providing Edge and AdjEntry sets.

## Symbol Visibilities
To align the behaviour of the Windows and UNIX versions, unix shared library builds now only make classes and functions
marked `OGDF_EXPORT` publicly visible (as is needed for Windows DLLs anyways). If a non-templated class or top-level
function you previously used (probably on Linux or MacOS) now became hidden,
please check the documentation of `OGDF_EXPORT` and then open a PR to make it visible.

## GraphIO

### SvgPrinter
When their height is equal to their width, `Shape::Triangle` and `Shape::InvTriangle` are now drawn as regular triangles in SVGs (and not as just isosceles ones).

### TikzWriter
`TikzWriter::isCoveredBy()` was removed in favor of `isPointCoveredByNode()` in `geometry.h`.

## RegisteredArrays
`RegisteredArray`, the underlying class for `NodeArray`, `EdgeArray` etc. now uses `std::vector` instead of `ogdf::Array` for data storage.
If the stored elements have a non-trivial move-constructor, it should be marked `noexcept`.
Otherwise, all elements will be [copied when the array grows](https://stackoverflow.com/a/28627764).

## NodeSet, FaceSet, and ClusterSet(Simple|Pure)
The template classes `NodeSet<bool SupportFastSizeQuery = true>` and `FaceSet<bool>` were converted to non-templated classes
using their default `SupportFastSizeQuery = true` versions, which now always keeps track of its size at a negligible overhead
(these sets do not support merging by splicing as simple double-linked lists do anyways).
Similarly, `ClusterSetPure` was removed in favor of `ClusterSet`.
`ClusterSetSimple` was removed in favor of `ClusterArray<bool>`.
All these `RegisteredSet`s now automatically remove members (nodes, faces, clusters,...) from their lists when they are deleted
from the corresponding registries (Graphs, CombEmbeddings, ClusterGraphs,...).

## Graph
The move constructor and assignment operators of `Graph` are now deleted, which especially means that `NodeArray<Graph>`, `EdgeArray<Graph>`, `FaceArray<Graph>`, etc. is no longer possible.
(Previously it compiled, but randomly broke at runtime when adding new nodes.)
Use `NodeArrayP<Graph>`, `EdgeArrayP<Graph>`, `FaceArrayP<Graph>`, etc. instead to wrap the `Graph`s in `std::unique_ptr`s and then make one pass over all entries to initialize the pointers.
See the documentation of `NodeArrayP` for usage as member variable on MSVC<=16.
For this reason, `SimDraw::getBasicGraph` now returns a `std::unique_ptr<GraphCopy>` instead of copying a `Graph` object on return.

## GraphCopy
`GraphCopy::createEmpty()` was deprecated in favor of `setOriginalGraph()`.
The same holds for `createEmpty()` of `GraphCopySimple` and `EdgeWeightedGraph`.

## (Hyper/Cluster)GraphObserver
`GraphObserver`s are now notified when their `Graph` is destructed through `GraphObserver::registrationChanged()`.

`ClusterGraphObserver`s are now notified of their `ClusterGraph` being cleared through `ClusterGraphObserver::clustersCleared()`.

`HypergraphObserver::init()` was deprecated in favor of `reregister()`.

`Observer`s and their `Observable`s now have deleted copy and move constructors and assignment operators.
Subclasses can instead explicitly declare their copy and move behaviour using the default constructors of `Observer` / `Observable`,
`Observer::getObservers()`, `Observer::clearObservers()` and `Observable::reregister()`.

Additionally, the `Observer(Observable*)` constructor (e.g. `GraphObserver(Graph*)`) is now deprecated as it would
trigger a `registrationChanged` callback before the construction of your subclass is done.
Instead, use the default `Observer()` constructor and explicitly call `reregister(...)` in the constructor of your subclass.

## Graph::insert
Multiple methods for inserting (parts of) a graph were merged into a single `Graph::insert` implementation.
This implementation is also used when copy-constructing or (in combination with clear) when copy-assigning.
It replaces the following different implementations, which were removed:

```c++
void Graph::construct(const Graph &G, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void Graph::assign(const Graph &G, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void Graph::insert(const Graph &G, NodeArray<node> &nodeMap);
void Graph::copy(const Graph &G);
void Graph::copy(const Graph &G, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void Graph::constructInitByNodes(const Graph &G, const List<node> &nodeList, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void Graph::constructInitByActiveNodes(const List<node> &nodeList, const NodeArray<bool> &activeNodes, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void Graph::constructInitByCC(const CCsInfo &info, int cc, NodeArray<node> &mapNode, EdgeArray<edge> &mapEdge);
void GraphCopySimple::init(const Graph &G);
void GraphCopySimple::initGC(const GraphCopySimple &GC, NodeArray<node> &vCopy, EdgeArray<edge> &eCopy);
void GraphCopy::init(const Graph &G);
void GraphCopy::initByCC(const CCsInfo &info, int cc, EdgeArray<edge> &eCopy);
void GraphCopy::initByNodes(const List<node> &origNodes, EdgeArray<edge> &eCopy);
void GraphCopy::initByActiveNodes(const List<node> &nodeList, const NodeArray<bool> &activeNodes, EdgeArray<edge> &eCopy);
void GraphCopy::initGC(const GraphCopy &GC, NodeArray<node> &vCopy, EdgeArray<edge> &eCopy);
// from <extended_graph_alg.h>:
void inducedSubGraph(const Graph &G, LISTITERATOR start, Graph &subGraph);
void inducedSubGraph(const Graph &G, LISTITERATOR start, Graph &subGraph, NodeArray<node> &nodeTableOrig2New);
void inducedSubGraph(const Graph &G, LISTITERATOR start, Graph &subGraph, NodeArray<node> &nodeTableOrig2New, EdgeArray<edge> &edgeTableOrig2New);
void inducedSubGraph(const Graph &G, LISTITERATOR start, GraphCopySimple &subGraph);
```

`GraphCopy::insert` (and `GraphCopySimple::insert`) will automatically update its mappings when inserting parts of the original Graph.
This can be disabled by using `setLinkCopiesOnInsert`.

## PoolMemoryAllocator

`PoolMemoryAllocator::defrag()` was renamed to `defragGlobal()` and (more importantly) now has a companion method `defragThread()`
that defragments the thread-local memory pool.

## Random Clusterings
The following methods were moved to `graph_generators/clustering.h` and renamed to more accurately reflect their functionality.
Note that `randomClusterPlanarGraph` did not actually generate cluster planar, but cluster-connected instances, and none of the methods used the `Graph` parameter.
```
randomClusterPlanarGraph(ClusterGraph& C, Graph& G, int cNum) -> randomCConnectedClustering(ClusterGraph& C, int cNum)
randomClusterGraph(ClusterGraph& C, Graph& G, int cNum) -> randomClustering(ClusterGraph& C, int cNum)
randomClusterGraph(ClusterGraph& C, const Graph& G, const node root, int moreInLeaves) -> randomClustering(ClusterGraph& C, const node root, int moreInLeaves)
```
