[OGDF](../../README.md) » [Release Notes](../relnotes.md) » Foxglove

# OGDF Foxglove (2025.10)

Released 2025-10-27.

This release contains many significant improvements to core functionality,
including a unified `Graph::insert` and `NodeArray`s/`EdgeArray`s that store
their data in `std::vector`s.
Moreover, we enhanced the build process by, e.g., properly enabling multi-config
builds and packaging via CPack.
However, as a result, this release also comes with several notable breaking
changes. For example, OGDF now requires C++17. Thus, we strongly suggest to read
the [porting guide](../porting/foxglove.md).

Note that to simplify access to the OGDF, we also plan to publish an
OGDF [docker image](https://hub.docker.com/r/ogdf/ogdf)
and to improve support of C++ package managers in the future.

Noteworthy changes:
 * build:
    * OGDF now requires C++17 features
    * improved compatibility with Windows MinGW, Apple silicon, and CGAL 6
    * fixes for linking user code with out-of-source OGDF builds
    * aligned symbol visibility on Windows and Unix
    * updated install declaration using GNUInstallDirs
    * enabled multi-config builds (i.e., parallel debug and release builds)
        * debug binaries and libs are now marked with `-debug`
    * enabled packaging via CPack
 * core functionality:
    * new unified and versatile `Graph::insert` implementation:
        * variants to insert subgraphs via iterators, lists, filter functions
        * options for notifying observers and keeping embedding or node/edge IDs
        * replaces `Graph::copy()`, `Graph::construct*()`, `Graph::assign()`,
          and `GraphCopy(Simple)::init*()`
    * new `RegisteredArray`:
        * as underlying class for `NodeArray`, `EdgeArray`, `FaceArray` etc.
        * now uses `std::vector` instead of `ogdf::Array` for data storage
        * new `NodeArrayP` etc. to wrap non-movable objects in `std::unique_ptr`
        * new `invertRegisteredArray` to transfer `XArray<Y>` into `YArray<X>`
    * new `RegisteredSet`:
        * as underlying class for `NodeSet`, `FaceSet`, `ClusterSet`
        * automatically removes members that are deleted from the resp. registry
        * always allows for efficient `size()`
        * new `EdgeSet` and `AdjEntrySet`
    * new `Observer` and `Observable`:
        * as underlying class for `GraphObserver`, `RegisteredSet`s etc.
        * new (`Graph`)`Observer::registrationChanged()`
        * new `ClusterGraphObserver::clustersCleared()`
    * `GraphCopy` and `GraphCopySimple`:
        * now with common superclass `GraphCopyBase`
        * new `GraphCopyBase::setOriginalGraph()` replacing `createEmpty()`
        * new `GraphCopySimple::copyEmbeddingToOriginal()`
        * new generic `copyEmbedding()`
    * `ClusterGraph`:
        * new `ClusterGraph::representsConnectedCombEmbedding()`
        * new `ClusterGraph::planarizeClusterBorderCrossings()`
        * new `ClusterGraph::adjAvailable()`
        * new `cluster->isDescendant()`
 * graph algorithms:
    * new algorithms for Min-Weight Perfect Matching (MWPM):
        * `MatchingModule`
        * `MatchingBlossom`
        * `MatchingBlossomV`
        * `MatchingImplementation` as the current best OGDF MWPM algorithm
    * new node coloring heuristics:
        * `NodeColoringModule`
        * `NodeColoringBergerRompel`
        * `NodeColoringBoppanaHalldorsson`
        * `NodeColoringHalldorsson`
        * `NodeColoringJohnson`
        * `NodeColoringRecursiveLargestFirst`
        * `NodeColoringSequential`
        * `NodeColoringSimple`
        * `NodeColoringWigderson`
    * new `PCTree` for planarity testing
    * new Synchronized Planarity functionality:
        * `SyncPlan` for modeling and solving SyncPlan instances
        * `ClusterPlanarityModule` for computing cluster-planar embeddings
        * `SyncPlanClusterPlanarityModule`
        * `randomSyncPlanInstance()` generator
        * `randomSEFEInstanceBy...()` generators
 * graph decomposition:
    * new `FourBlockTree` for construction and traversal of 4-block trees
    * new parameter for `connectedComponents()` to get comp representatives
    * new `CCsInfo::nodes(int cc)` and `CCsInfo::edges(int cc)`
    * new access methods for member vars of `BCTree` and `DynamicSPQRForest`
    * new `operator<<` for nodes of `BCTree` and `DynamicSPQRForest`
 * graph generators:
    * new `randomProperMaximalLevelPlaneGraph()`
    * new `pruneEdges()` to enforce a maximum edge number
    * new graph operations `join()`, `intersection()`, and `complement()`
    * reworked cluster generators:
        * old `randomClusterPlanarGraph`  renamed `randomCConnectedClustering()`
        * old `randomClusterGraph`  renamed `randomClustering()`
        * new `randomClusterPlanarGraph`
        * new `randomPlanarClustering`
 * layouts:
    * `TreeLayout` now works on non-arborescence forests
    * `SimpleCCPacker` now correctly preserves edge bends
    * `ComponentSplitterLayout` now correctly preserves edge bends
 * memory allocation:
    * `PoolMemoryAllocator::defrag()` renamed `defragGlobal()`
    * new `PoolMemoryAllocator::defragThread()`
    * new `PoolMemoryAllocator::get{Global,Thread}FreeListSizes()`
    * new `OGDFAllocator` for use with containers of the C++ standard lib
 * miscellaneous:
    * `HiddenEdgeSet` now has `begin()`, `end()`, and `empty()`
    * new `Graph::sort(node v, ITER begin, ITER end)` to adapt rotation systems
    * new indent parameter for `Logger::lout()`
        * new `Logger::setIndent()`, `Logger::indent()`, `Logger::dedent()`
    * `SvgPrinter` now correctly connects arrow heads to the edge's target
 * OGDF development:
    * documentation now powered by Netlify
    * new `CODE_OF_CONDUCT.md`
    * new utility `OGDF_`-macros for declaring c'tors and assignment operators
    * new `OGDF_IF_DBG` macro for single-line statements in debug mode
    * new `indent_comments.py` for formatting comments
    * new `style/test_all.sh` for unified style check
    * usage of [include-what-you-use](https://include-what-you-use.org/)
    * new `make_release.sh` for easier release management

This release contains (big and small) contributions by Antoine Lambert, Dominik
Potulski, Felix F Xu, Gregor Diatzko, Jan-Niklas Buckow, Joshua Sangmeister,
Lily Wang, Max Ilsen, Simon Dominik "Niko" Fink, and Sven Strickroth.
Thanks a lot to all contributors!
