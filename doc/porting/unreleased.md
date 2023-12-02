[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Unreleased

# Porting from Elderberry to current unreleased version

## GraphIO

### SvgPrinter
When their height is equal to their width, `Shape::Triangle` and `Shape::InvTriangle` are now drawn as regular triangles in SVGs (and not as just isosceles ones).

### TikzWriter
`TikzWriter::isCoveredBy()` was removed in favor of `isPointCoveredByNode()` in `geometry.h`.

## RegisteredArrays and Graph::insert

`RegisteredArray` now uses `std::vector` instead of `ogdf::Array` for data storage.
If the store elements have a non-trivial move-constructor, it should be marked `noexcept`.
Otherwise, all elements will be [copied when the array grows](https://stackoverflow.com/a/28627764).
Removed `ClusterSetSimple` in favour of `ClusterArray<bool>` and `ClusterSetPure` in favour of `ClusterSet<false>` (which does not keep track of its size).


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

`GraphObserver`s are now notified when their `Graph` is destructed through `GraphObserver::registrationChanged()`.
