[OGDF](/README.md) » [Porting Guide](/doc/porting.md) » Snapshot

# Porting from Baobab to current snapshot

## General

### C++11

OGDF now requires C++11 features.
Make sure you do not use another standard when compiling your user programs.

### Compiling for debug mode

If you want to compile for debug mode, you do not need to define
`OGDF_DEBUG` any longer. This is automatically defined using
CMake.

### Header files

It might be necessary that you need to include additional header files
for code that worked before. For example, `#include <NodeArray.h>`
does not ensure that you also have an `EdgeArray`.

### Global namespace

In former versions of the OGDF some symbols were added to the global namespace.
This includes but isn't limited to `ifstream`, `ofstream`, `min`, `max`, and `numeric_limits`.
Thus, you may be required to explicitly use the OGDF namesapce in some places.

## Iteration Macros

We have removed most of the iteration macros since C++11 offers range-based for loops.
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
