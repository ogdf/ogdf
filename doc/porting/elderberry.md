[OGDF](../../README.md) » [Developer's Guide](../dev-guide.md) » [Porting Guide](../porting.md) » Elderberry

# Porting from Dogwood to Elderberry

## Compilation: `OGDF_LEAK_CHECK` affecting COIN

`OGDF_LEAK_CHECK` enables the address sanitizer for COIN as well.
This prevents false positives resulting from different parts of an application being built with and without the address sanitizer.
See [here](https://github.com/google/sanitizers/wiki/AddressSanitizerContainerOverflow#false-positives).

## GraphObserver

`GraphObserver::reInit()` was removed as it was never actually called.

## MinCut

The class `MinCut` was renamed to `MinimumCutStoerWagner` and is now a template,
e.g., you can use `MinimumCutStoerWagner<int>` when you have integer weights.
Its methods were also renamed and now return const-references (where appropriate)
instead of copying data.
The returned lists are now ArrayBuffers.

The code
```c++
MinCut mincut(graph, weights);
double value = mincut.minimumCut(); // computation here
double value_again = mincut.minCutValue();
List<node> nodes;
mincut.partition(nodes);
std::cout << nodes << std::endl;
List<edge> edges;
mincut.cutEdges(edges, graph);
std::cout << edges << std::endl;
```
transforms to
```c++
MinimumCutStoerWagner<int> mincut;
double value = mincut.call(graph, weights); // computation here
double value_again = mincut.value();
std::cout << mincut.nodes() << std::endl;
std::cout << mincut.edges() << std::endl;
```

## RadialTreeLayout

`RadialTreeLayout::leaves()` and `RadialTreeLayout::connectedComponentDistance()` have been removed.
The latter was not properly implemented and had no effect on the layout.
